"""Load BiGG JSON or metabolic SBML with fixed stoichiometry."""
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def sbml_metadata(element):
    """Keep direct annotations and qualifiers, without borrowing child annotations."""
    annotations = []
    notes = ''
    for child in element:
        local = child.tag.rsplit('}', 1)[-1]
        if local == 'notes':
            notes = ' '.join(' '.join(child.itertext()).split())
        if local != 'annotation':
            continue
        for description in child.iter():
            if description.tag != '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description':
                continue
            if description.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about') != '#' + element.get('metaid', ''):
                continue
            for qualifier in description:
                if not qualifier.tag.startswith(('{http://biomodels.net/biology-qualifiers/}',
                                                 '{http://biomodels.net/model-qualifiers/}')):
                    continue
                for resource in qualifier.iter():
                    uri = resource.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource')
                    if uri:
                        annotations.append(dict(qualifier=qualifier.tag, uri=uri))
    return dict(annotations=annotations, notes=notes)


def load_model(path):
    path = Path(path)
    if path.suffix.lower() == '.json':
        return json.loads(path.read_text(encoding='utf-8'))
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise ValueError(f'Invalid SBML XML: {exc}') from exc
    if root.tag.split('}')[-1] != 'sbml':
        raise ValueError('Expected an SBML document')
    namespace = root.tag.rsplit('}', 1)[0] + '}' if '}' in root.tag else ''
    def children(element, name):
        return element.findall('/'.join(namespace + p for p in name.split('/')))
    models = children(root, 'model')
    if len(models) != 1:
        raise ValueError('SBML must contain one model')
    model = models[0]
    if any(children(model, name) for name in ('listOfRules', 'listOfEvents', 'listOfInitialAssignments')):
        raise ValueError('Dynamic SBML rules, events and initial assignments are unsupported')
    def decode(identifier, prefix):
        if identifier.startswith(prefix):
            identifier = identifier[len(prefix):]
        return re.sub(r'__(\d+)__', lambda m: chr(int(m[1])), identifier)
    species = children(model, 'listOfSpecies/species')
    id_map = {s.attrib['id']: decode(s.attrib['id'], 'M_') for s in species}
    if len(id_map) != len(species) or len(set(id_map.values())) != len(species):
        raise ValueError('Duplicate metabolite IDs after SBML ID normalization')
    metabolites = [dict(id=id_map[s.attrib['id']], name=s.get('name', ''),
                        compartment=s.get('compartment', ''), **sbml_metadata(s),
                        **{k.rsplit('}', 1)[-1]: v for k, v in s.attrib.items()
                           if k.endswith(('}chemicalFormula', '}charge'))}) for s in species]
    def fbc_tag(element, name):
        return '/fbc/' in element.tag and element.tag.rsplit('}', 1)[-1] == name

    def fbc_attribute(element, name):
        return next((v for k, v in element.attrib.items()
                     if '/fbc/' in k and k.rsplit('}', 1)[-1] == name), element.get(name, ''))

    genes = [dict(id=decode(fbc_attribute(g, 'id'), 'G_'), name=fbc_attribute(g, 'name'),
                  label=fbc_attribute(g, 'label'), **sbml_metadata(g))
             for g in model.iter() if fbc_tag(g, 'geneProduct')]

    def gene_rule(node):
        tag = node.tag.rsplit('}', 1)[-1]
        if tag == 'geneProductRef':
            return decode(fbc_attribute(node, 'geneProduct'), 'G_')
        parts = [gene_rule(child) for child in node]
        if tag in ('and', 'or'):
            return '(' + (' ' + tag + ' ').join(parts) + ')'
        return ' '.join(parts)
    parameters = {p.attrib['id']: float(p.attrib['value'])
                  for p in children(model, 'listOfParameters/parameter') if 'value' in p.attrib}
    reactions = []
    for reaction in children(model, 'listOfReactions/reaction'):
        rid = decode(reaction.attrib['id'], 'R_')
        coefficients = {}
        for group, sign in (('listOfReactants', -1), ('listOfProducts', 1)):
            for ref in children(reaction, group + '/speciesReference'):
                if children(ref, 'stoichiometryMath') or ref.get('constant') == 'false':
                    raise ValueError(f'Variable stoichiometry is unsupported: {rid}')
                raw = ref.attrib['species']
                if raw not in id_map:
                    raise ValueError(f'Unknown species {raw} in reaction {rid}')
                value = float(ref.get('stoichiometry', '1'))
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f'Invalid stoichiometry in reaction {rid}')
                mid = id_map[raw]
                coefficients[mid] = coefficients.get(mid, 0) + sign * value
        # COBRA Level 2 stores bounds as kinetic-law parameters.
        legacy = {p.attrib['id']: float(p.attrib['value']) for p in
                  children(reaction, 'kineticLaw/listOfParameters/parameter') if 'value' in p.attrib}
        bounds = [legacy.get('LOWER_BOUND', -1000 if reaction.get('reversible', 'true') in ('true', '1') else 0),
                  legacy.get('UPPER_BOUND', 1000)]
        for i, name in enumerate(('lowerFluxBound', 'upperFluxBound')):
            refs = [v for k, v in reaction.attrib.items()
                    if k.startswith('{') and '/fbc/' in k and k.split('}')[-1] == name]
            if refs:
                if refs[0] not in parameters:
                    raise ValueError(f'Unknown flux-bound parameter {refs[0]} in {rid}')
                bounds[i] = parameters[refs[0]]
        if any(math.isnan(b) for b in bounds) or bounds[0] > bounds[1]:
            raise ValueError(f'Invalid flux bounds in reaction {rid}')
        associations = [c for c in reaction if fbc_tag(c, 'geneProductAssociation')]
        reactions.append(dict(id=rid, name=reaction.get('name', ''), **sbml_metadata(reaction),
                              gene_reaction_rule=gene_rule(associations[0]) if associations else '',
                              metabolites={m: c for m, c in coefficients.items() if c},
                              lower_bound=bounds[0], upper_bound=bounds[1]))
    if not metabolites or not reactions or len({r['id'] for r in reactions}) != len(reactions):
        raise ValueError('Model must contain metabolites and reactions with unique IDs')
    return dict(id=model.get('id', ''), name=model.get('name', ''), **sbml_metadata(model),
                metabolites=metabolites, reactions=reactions, genes=genes)


def model_network(model):
    """Use the same undirected co-participation graph as the CSV generator."""
    from stoichiometric_network import stoichiometric_matrix, adjacency_from_stoichiometry
    matrix, ids, _ = stoichiometric_matrix(model)
    adjacency = adjacency_from_stoichiometry(matrix)
    return ids, [(i, j) for i in range(len(ids)) for j in range(i + 1, len(ids))
                 if adjacency[i][j]]
