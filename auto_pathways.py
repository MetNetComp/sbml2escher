"""Evidence-based pathway discovery. External data never replaces model reactions.

KEGG KGML describes representative substrate/product pairs, not complete equations.
Pair matches therefore remain candidates, even when gene evidence supports them.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, urlparse, unquote
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

KEGG = 'https://rest.kegg.jp/'
NCBI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
CHEBI = 'https://www.ebi.ac.uk/chebi/backend/api/public/compound/'
REACTOME = 'https://reactome.org/ContentService/data/'
COFACTORS = {'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006',
             'C00007', 'C00008', 'C00009', 'C00011', 'C00013', 'C00016',
             'C00010', 'C00020', 'C00080', 'C01352'}
CENTRAL = {'00010', '00020', '00030', '00190', '00195'}
DEFAULT_WEIGHTS = dict(coverage=0.25, continuity=0.25, evidence=0.20, central=0.30)


class DiscoveryError(ValueError):
    pass


class DatabaseClient:
    """Content-addressed responses, explicit offline misses, bounded retries."""
    def __init__(self, cache, offline=False, refresh=False, timeout=25):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline, self.refresh, self.timeout = offline, refresh, timeout
        self.requests, self.last_request = {}, 0

    def get(self, url):
        key = hashlib.sha256(url.encode()).hexdigest()
        path = self.cache / (key + '.json')
        if path.exists() and (self.offline or not self.refresh):
            record = json.loads(path.read_text())
            if record['url'] != url or hashlib.sha256(record['body'].encode()).hexdigest() != record['sha256']:
                raise DiscoveryError('Corrupt database cache: ' + str(path))
        else:
            if self.offline:
                raise DiscoveryError('Offline cache miss: ' + url)
            for attempt in range(3):
                time.sleep(max(0, 0.4 - (time.monotonic() - self.last_request)))
                self.last_request = time.monotonic()
                try:
                    with urlopen(Request(url, headers={'User-Agent': 'sbml2escher/auto-pathways',
                                                       'Accept': '*/*'}), timeout=self.timeout) as response:
                        body = response.read().decode('utf-8')
                    break
                except (HTTPError, URLError, TimeoutError) as exc:
                    if isinstance(exc, HTTPError) and exc.code not in (429, 500, 502, 503, 504):
                        raise DiscoveryError(f'Database request failed: {url}: {exc}') from exc
                    if attempt == 2:
                        raise DiscoveryError(f'Database request failed: {url}: {exc}') from exc
                    time.sleep(2 ** attempt)
            record = dict(url=url, retrieved_at=datetime.now(timezone.utc).isoformat(),
                          sha256=hashlib.sha256(body.encode()).hexdigest(), body=body)
            temporary = path.with_suffix('.tmp')
            temporary.write_text(json.dumps(record, ensure_ascii=False))
            temporary.replace(path)
        self.requests[url] = {k: v for k, v in record.items() if k != 'body'}
        return record['body']

    def json(self, url):
        try:
            return json.loads(self.get(url))
        except json.JSONDecodeError as exc:
            raise DiscoveryError('Invalid JSON response: ' + url) from exc


def annotation_ids(entity, database, qualifiers=('is',)):
    """Identifiers.org URI and BiGG JSON annotations, restricted by qualifier."""
    result = set()
    for item in entity.get('annotations', []):
        if item['qualifier'].rsplit('}', 1)[-1] not in qualifiers:
            continue
        uri = unquote(item['uri'])
        if uri.startswith('urn:miriam:'):
            namespace, _, value = uri[len('urn:miriam:'):].partition(':')
        else:
            parsed = urlparse(uri)
            if parsed.hostname not in ('identifiers.org', 'www.identifiers.org'):
                continue
            part = parsed.path.strip('/')
            if '/' in part:
                namespace, value = part.split('/', 1)
            else:
                namespace, _, value = part.partition(':')
        if namespace == database and value:
            result.add(value)
    values = entity.get('annotation', {}).get(database, [])
    if isinstance(values, str):
        values = [values]
    result.update(values)
    return sorted(result)


def taxonomy(client, taxid):
    data = client.json(NCBI + 'esummary.fcgi?' + urlencode(dict(db='taxonomy', id=taxid, retmode='json')))
    record = data.get('result', {}).get(str(taxid), {})
    if not record.get('scientificname'):
        raise DiscoveryError('Taxonomy ID was not resolved: ' + str(taxid))
    return dict(taxid=str(taxid), scientific_name=record['scientificname'],
                rank=record.get('rank'), source='NCBI Taxonomy')


def infer_organism(model, client, override=None, kegg_override=None):
    evidence, warnings = [], []
    explicit = annotation_ids(model, 'taxonomy', ('is', 'hasTaxon'))
    explicit += annotation_ids(model, 'ncbitaxon', ('is', 'hasTaxon'))
    explicit = sorted(set(explicit))
    if override:
        if override.isdigit():
            choice = taxonomy(client, override)
        else:
            term = override + '[Scientific Name]'
            found = client.json(NCBI + 'esearch.fcgi?' + urlencode(dict(db='taxonomy', term=term, retmode='json')))
            ids = found.get('esearchresult', {}).get('idlist', [])
            if len(ids) != 1:
                raise DiscoveryError('Organism name is ambiguous or unknown; specify a Taxonomy ID')
            choice = taxonomy(client, ids[0])
        status = 'user_specified'
        evidence.append(dict(source='command_line', value=override))
        if explicit and choice['taxid'] not in explicit:
            warnings.append('User organism differs from model taxonomy annotation')
    elif len(explicit) == 1:
        choice = taxonomy(client, explicit[0])
        status = 'explicit'
        evidence.append(dict(source='model_annotation', taxid=explicit[0]))
    elif explicit:
        return dict(status='ambiguous', candidates=explicit, evidence=evidence,
                    warnings=['Multiple model taxa; no single organism selected'])
    else:
        # Cross-check several referenced genes; shared metabolites do not identify taxa.
        ids = sorted({v for gene in model.get('genes', []) for db in ('ncbigene', 'ncbi.gene')
                      for v in annotation_ids(gene, db) if v.isdigit()})[:5]
        organisms = {}
        if ids:
            try:
                data = client.json(NCBI + 'esummary.fcgi?' + urlencode(dict(db='gene', id=','.join(ids), retmode='json')))
                for gid in ids:
                    organism = data.get('result', {}).get(gid, {}).get('organism', {})
                    if organism.get('taxid'):
                        tid = str(organism['taxid'])
                        organisms[tid] = organism
                        evidence.append(dict(source='NCBI Gene', gene_id=gid, taxid=tid,
                                             scientific_name=organism.get('scientificname')))
            except DiscoveryError as exc:
                warnings.append(str(exc))
        if len(organisms) == 1 and len(evidence) >= 2:
            choice = taxonomy(client, next(iter(organisms)))
            status = 'inferred_from_genes'
        elif len(organisms) > 1:
            return dict(status='ambiguous', candidates=sorted(organisms), evidence=evidence,
                        warnings=warnings + ['Referenced genes belong to multiple taxa'])
        else:
            text = ' '.join([model.get('id', ''), model.get('name', ''), model.get('notes', '')])
            names = {'e_coli': 'Escherichia coli', 'Escherichia coli': 'Escherichia coli',
                     'Homo sapiens': 'Homo sapiens', 'Saccharomyces cerevisiae': 'Saccharomyces cerevisiae'}
            matches = sorted({name for token, name in names.items() if token.lower() in text.lower()})
            matches += re.findall(r'(?:organism|species|生物種)\s*[:=]\s*([A-Z][a-z]+ [a-z]+)\b', text)
            matches = sorted(set(matches))
            if len(matches) != 1:
                return dict(status='unresolved', candidates=matches, evidence=evidence, warnings=warnings)
            found = client.json(NCBI + 'esearch.fcgi?' + urlencode(dict(
                db='taxonomy', term=matches[0] + '[Scientific Name]', retmode='json')))
            ids = found.get('esearchresult', {}).get('idlist', [])
            if len(ids) != 1:
                return dict(status='unresolved', candidates=ids, evidence=evidence, warnings=warnings)
            choice = taxonomy(client, ids[0])
            status = 'inferred_from_model_text'
            evidence.append(dict(source='model_text', scientific_name=matches[0]))
    result = dict(choice, status=status, evidence=evidence, warnings=warnings)
    genomes = []
    for line in client.get(KEGG + 'list/genome/' + choice['taxid']).splitlines():
        fields = line.split('\t')
        if len(fields) >= 2:
            code = fields[1].split(';')[0].strip()
            genomes.append(dict(genome_id=fields[0].removeprefix('gn:'), code=code,
                                description=fields[-1]))
    result['kegg_candidates'] = genomes
    if kegg_override:
        if not re.fullmatch(r'[a-z][a-z0-9]{2,4}', kegg_override):
            raise DiscoveryError('Invalid KEGG organism code')
        if kegg_override not in {g['code'] for g in genomes}:
            raise DiscoveryError('KEGG organism code does not belong to the resolved taxon')
        result['kegg_organism'] = kegg_override
    elif len(genomes) == 1:
        result['kegg_organism'] = genomes[0]['code']
    else:
        result['kegg_organism'] = None
        result['warnings'].append('No unique KEGG genome; using reference pathways without a strain claim')
    return result


def compound_mapping(model):
    records = []
    for m in model['metabolites']:
        ids = [v for v in annotation_ids(m, 'kegg.compound') if re.fullmatch(r'C\d{5}', v)]
        records.append(dict(model_id=m['id'], name=m.get('name', ''), compartment=m.get('compartment', ''),
                            kegg_compounds=ids, chebi_ids=annotation_ids(m, 'chebi'),
                            status='unique' if len(ids) == 1 else 'ambiguous' if ids else 'unmapped',
                            source='SBML/JSON annotation',
                            formula=m.get('chemicalFormula', m.get('formula')), charge=m.get('charge')))
    return records


def pairs(text):
    return [line.split('\t')[:2] for line in text.splitlines() if '\t' in line]


def batches(items, size=10):
    items = sorted(set(items))
    for i in range(0, len(items), size):
        yield items[i:i + size]


def parse_kgml(text, organism=None):
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise DiscoveryError('Invalid KEGG KGML') from exc
    if root.tag != 'pathway':
        raise DiscoveryError('KEGG response is not a KGML pathway')
    entries = {e.get('id'): e for e in root.findall('entry')}
    records = []
    for reaction in root.findall('reaction'):
        entry = entries.get(reaction.get('id'))
        if organism and (entry is None or entry.get('type') != 'gene'):
            continue  # KGML also contains ortholog boxes absent in this organism.
        genes = entry.get('name', '').split() if entry is not None else []
        for kid in reaction.get('name', '').split():
            if not re.fullmatch(r'rn:R\d{5}', kid):
                continue
            records.append(dict(kegg_reaction=kid[3:], genes=genes,
                                substrates=sorted({c for s in reaction.findall('substrate')
                                                   for c in s.get('name', '').split() if c.startswith('cpd:')}),
                                products=sorted({c for s in reaction.findall('product')
                                                 for c in s.get('name', '').split() if c.startswith('cpd:')}),
                                reversible=reaction.get('type') == 'reversible'))
    return records


def automatic_currency(model, mappings):
    """DB identities fix name-based misses (e.g. 'ATP C10H12N5O13P3')."""
    names = {'h', 'h+', 'h2o', 'water', 'o2', 'oxygen', 'co2', 'atp', 'adp', 'amp',
             'pi', 'ppi', 'nad', 'nadh', 'nadp', 'nadph', 'fad', 'fadh2'}
    return {m['model_id'] for m in mappings if
            (m['kegg_compounds'] and set(m['kegg_compounds']) <= COFACTORS) or
            re.sub(r'\[[^]]+\]$|_[a-z][0-9]?$', '', m['name']).lower() in names}


def match_pathway(model, db_reactions, mappings, currency, organism):
    compounds = {m['model_id']: set(m['kegg_compounds']) for m in mappings}
    genes = {g['id']: g for g in model.get('genes', [])}
    matches = []
    for reaction in model['reactions']:
        if reaction['lower_bound'] == reaction['upper_bound'] == 0:
            continue
        annotated = set(annotation_ids(reaction, 'kegg.reaction'))
        rule_ids = set(re.findall(r'[A-Za-z0-9_.:-]+', reaction.get('gene_reaction_rule', ''))) - {'and', 'or'}
        gene_ids = set()
        for gid in rule_ids:
            g = genes.get(gid, {})
            gene_ids.update(annotation_ids(g, 'kegg.genes'))
            if organism:
                gene_ids.update(organism + ':' + v for v in annotation_ids(g, 'refseq_locus_tag'))
        options = []
        for db in db_reactions:
            direct = db['kegg_reaction'] in annotated
            gene_match = sorted(gene_ids & set(db['genes']))
            orientations = [(db['substrates'], db['products'])]
            if db['reversible']:
                orientations.append((db['products'], db['substrates']))
            for direction, allowed in [('forward', reaction['upper_bound'] > 0),
                                       ('reverse', reaction['lower_bound'] < 0)]:
                if not allowed:
                    continue
                for a, ca in reaction['metabolites'].items():
                    if a in currency or not compounds.get(a) or (ca < 0) != (direction == 'forward'):
                        continue
                    for b, cb in reaction['metabolites'].items():
                        if b in currency or not compounds.get(b) or (cb > 0) != (direction == 'forward'):
                            continue
                        compatible = [(compounds[a] & {v[4:] for v in aa}, compounds[b] & {v[4:] for v in bb})
                                      for aa, bb in orientations]
                        compatible = [(aa, bb) for aa, bb in compatible if len(aa) == len(bb) == 1]
                        if len(compatible) != 1:
                            continue
                        aa, bb = compatible[0]
                        # A unique intersection resolves the incidence, not the global identity.
                        confidence = 1.0 if direct else 0.8 if gene_match else 0.35
                        reason = 'reaction_id_and_participant_pair' if direct else 'gene_and_participant_pair' if gene_match else 'participant_pair_only'
                        options.append(dict(reaction_id=reaction['id'], input=a, output=b, direction=direction,
                                            kegg_reaction=db['kegg_reaction'], confidence=confidence,
                                            evidence=reason, matching_genes=gene_match,
                                            input_kegg_compound=next(iter(aa)), output_kegg_compound=next(iter(bb)),
                                            context_resolved_compounds=[mid for mid in (a, b) if len(compounds[mid]) > 1],
                                            stoichiometry_identity_verified=False))
        if options:
            best = sorted(options, key=lambda o: (-o['confidence'], o['kegg_reaction'], o['input'], o['output'], o['direction']))[0]
            best['alternative_pairs'] = len(options) - 1
            matches.append(best)
    return matches


def continuity(matches):
    """Largest connected reaction component; separate from database coverage."""
    adjacency = defaultdict(set)
    for m in matches:
        a, b, r = 'met:' + m['input'], 'met:' + m['output'], 'rxn:' + m['reaction_id']
        adjacency[a].add(r); adjacency[r].add(a)
        adjacency[b].add(r); adjacency[r].add(b)
    seen, components = set(), []
    for key in sorted(adjacency):
        if key in seen:
            continue
        todo, component = [key], set()
        while todo:
            node = todo.pop()
            if node in component:
                continue
            component.add(node); seen.add(node); todo.extend(adjacency[node] - component)
        components.append(sorted(n[4:] for n in component if n.startswith('rxn:')))
    return dict(largest_component_fraction=max(map(len, components), default=0) / max(len(matches), 1),
                reaction_components=components,
                note='Structural connectivity with allowed model directions; flux feasibility is not tested')


def backbone_layout(groups):
    """Condense directed SCCs; arrange branches in layers and cycles on circles."""
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    vertices, edges = set(), []
    for g in groups:
        for s in g['steps']:
            a, b, r = 'met:' + s['input'], 'met:' + s['output'], 'rxn:' + s['reaction_id']
            vertices.update((a, b, r)); edges.extend(((a, r), (r, b)))
    keys = sorted(vertices)
    index = {key: i for i, key in enumerate(keys)}
    matrix = coo_matrix((np.ones(len(edges)), ([index[a] for a, b in edges],
                                              [index[b] for a, b in edges])), shape=(len(keys), len(keys))).tocsr()
    count, labels = connected_components(matrix, directed=True, connection='strong')
    components = [sorted(key for key in keys if labels[index[key]] == i) for i in range(count)]
    outgoing, incoming = defaultdict(set), defaultdict(set)
    for a, b in edges:
        i, j = int(labels[index[a]]), int(labels[index[b]])
        if i != j:
            outgoing[i].add(j); incoming[j].add(i)
    rank, pending = {}, set(range(count))
    while pending:
        ready = sorted((i for i in pending if incoming[i] <= rank.keys()), key=lambda i: components[i])
        for i in ready:
            rank[i] = max((rank[j] + 1 for j in incoming[i]), default=0)
            pending.remove(i)
    layers = defaultdict(list)
    for i, r in rank.items():
        layers[r].append(i)
    coordinates, x = {}, 0.0
    for r in sorted(layers):
        layer = sorted(layers[r], key=lambda i: components[i])
        width = max(18, max(8 * len(components[i]) for i in layer))
        y = 0.0
        for i in layer:
            nodes = components[i]
            radius = max(12, 3 * len(nodes)) if len(nodes) > 1 else 0
            for k, node in enumerate(nodes):
                theta = 2 * math.pi * k / len(nodes)
                coordinates[node] = [x + width / 2 + radius * math.cos(theta), y + radius * math.sin(theta)]
            y -= 2 * radius + 24
        x += width + 24
    for g in groups:
        g['metabolite_positions'] = {mid: coordinates['met:' + mid] for s in g['steps']
                                    for mid in (s['input'], s['output'])}
        for s in g['steps']:
            s['position'] = coordinates['rxn:' + s['reaction_id']]
        p = min(g['metabolite_positions'].values(), key=lambda p: (p[0], -p[1]))
        g['label_position'] = [p[0], p[1] + 10]
    return [dict(nodes=c, kind='directed_structural_cycle_component') for c in components if len(c) > 1]


def important_compounds(model, mappings, selected, cycles, currency):
    membership = defaultdict(set)
    for candidate in selected:
        for step in candidate['selected_steps']:
            membership[step['input']].add(candidate['id'])
            membership[step['output']].add(candidate['id'])
    cyclic = {key[4:] for c in cycles for key in c['nodes'] if key.startswith('met:')}
    reactions = defaultdict(list)
    for r in model['reactions']:
        for mid in r['metabolites']:
            reactions[mid].append(r)
    records = []
    for mapping in mappings:
        mid = mapping['model_id']
        reasons = []
        if membership[mid]:
            reasons.append('selected_pathway_intermediate')
        if len(membership[mid]) > 1:
            reasons.append('pathway_connector')
        if mid in cyclic:
            reasons.append('directed_cycle_component')
        if mid in currency:
            reasons.append('energy_redox_or_auxiliary_compound')
        if any(len(r['metabolites']) == 1 for r in reactions[mid]):
            reasons.append('exchange_or_boundary_participant')
        if any('biomass' in (r['id'] + ' ' + r.get('name', '')).lower() for r in reactions[mid]):
            reasons.append('biomass_participant')
        if reasons:
            records.append(dict(mapping, reasons=reasons, pathways=sorted(membership[mid]),
                                pathway_importance=2 * len(membership[mid]) + (3 if mid in cyclic else 0),
                                auxiliary_importance=int(mid in currency), chebi_roles=[]))
    return sorted(records, key=lambda r: (-r['pathway_importance'], r['model_id']))


def enrich_chebi(records, client, warnings, limit=20):
    """Bounded enrichment; ambiguous chemical identities remain explicit candidates."""
    queried = {}
    ordered = sorted(records, key=lambda r: (-r['auxiliary_importance'], -r['pathway_importance'], r['model_id']))
    identifiers = [sorted(set(record['chebi_ids'])) for record in ordered]
    # Share the request budget across compounds before querying second alternatives.
    for depth in range(max(map(len, identifiers), default=0)):
        for ids in identifiers:
            if depth >= len(ids):
                continue
            cid = ids[depth]
            if not re.fullmatch(r'CHEBI:\d+', cid) or cid in queried or len(queried) >= limit:
                continue
            try:
                queried[cid] = client.json(CHEBI + quote(cid, safe=':') + '/')
            except DiscoveryError as exc:
                warnings.append(str(exc)); queried[cid] = None
    for record, ids in zip(ordered, identifiers):
        record['chebi_identity_status'] = 'unique' if len(ids) == 1 else 'ambiguous' if ids else 'unmapped'
        record['chebi_enrichment'] = 'not_requested'
        for cid in ids:
            data = queried.get(cid)
            if data:
                record['chebi_enrichment'] = 'retrieved'
                record['chebi_roles'].append(dict(annotation_id=cid, resolved_id=data.get('chebi_accession'),
                                                  chemical_name=data.get('name'),
                                                  roles=data.get('roles_classification', []),
                                                  modified_on=data.get('modified_on'),
                                                  annotation_unambiguous=len(ids) == 1))


def reactome_candidates(model, mappings, taxon, client, warnings, limit=20):
    """Supplementary species-specific membership; never infer reaction identity from it."""
    result = dict(source='Reactome', taxid=taxon['taxid'], candidates=[], mappings=[],
                  interpretation='Identifier membership is supporting evidence, not reaction identity or pathway activity',
                  query_limit=limit)
    try:
        species = client.json(REACTOME + 'species/all')
        supported = [s for s in species if str(s.get('taxId')) == taxon['taxid']]
        if not supported:
            result['status'] = 'unsupported_taxon'
            warnings.append('Reactome does not list the resolved taxon; no human projection or strain substitution performed')
            return result
        result['species_name'] = supported[0].get('displayName')
        identifiers = []
        for gene in model.get('genes', []):
            for kid in annotation_ids(gene, 'uniprot'):
                identifiers.append(('UniProt', kid, gene['id'], 'gene'))
        for m in mappings:
            # Do not use broad or ambiguous chemical references as membership proof.
            if len(m['chebi_ids']) == 1 and re.fullmatch(r'CHEBI:\d+', m['chebi_ids'][0]):
                identifiers.append(('ChEBI', m['chebi_ids'][0].split(':')[1], m['model_id'], 'metabolite'))
        by_identifier = defaultdict(list)
        for db, kid, mid, kind in identifiers:
            by_identifier[db, kid].append(dict(model_id=mid, kind=kind))
        candidates = {}
        for (db, kid), owners in sorted(by_identifier.items())[:limit]:
            url = REACTOME + 'mapping/' + db + '/' + quote(kid, safe='') + '/pathways?' + urlencode(dict(species=taxon['taxid']))
            pathways = client.json(url)
            if not isinstance(pathways, list):
                raise DiscoveryError('Unexpected Reactome mapping response: ' + url)
            ids = []
            for pathway in pathways:
                pid = pathway.get('stId')
                if not pid:
                    continue
                # The API species filter is mandatory; additionally reject a named mismatch.
                names = pathway.get('speciesName', [])
                if isinstance(names, str):
                    names = [names]
                if names and result['species_name'] not in names:
                    continue
                ids.append(pid)
                candidate = candidates.setdefault(pid, dict(id=pid, name=pathway.get('displayName'),
                    source='Reactome', source_url='https://reactome.org/content/detail/' + pid,
                    taxid=taxon['taxid'], evidence=[], inferred=pathway.get('isInferred'),
                    selected_for_backbone=False))
                candidate['evidence'].append(dict(database=db, identifier=kid, model_entities=owners))
            result['mappings'].append(dict(database=db, identifier=kid, model_entities=owners, pathway_ids=sorted(ids)))
        result['candidates'] = sorted(candidates.values(), key=lambda c: (-len(c['evidence']), c['id']))
        result['status'] = 'complete'
    except DiscoveryError as exc:
        warnings.append(str(exc))
        result['status'] = 'request_failed'
    return result


def discover(model, output, client, organism=None, kegg_organism=None, databases=('kegg',),
             max_pathways=5, candidate_limit=12, min_reactions=2, chebi_limit=20,
             weights=None, progress=print, reactome_limit=20):
    """Write auditable candidates and an executable, automatically placed backbone."""
    if 'kegg' not in databases or set(databases) - {'kegg', 'chebi', 'reactome'}:
        raise DiscoveryError('Supported databases: kegg, optional chebi and reactome')
    if not 1 <= max_pathways <= candidate_limit <= 100 or min_reactions < 1 or chebi_limit < 0 or reactome_limit < 0:
        raise DiscoveryError('Require 1 <= max-pathways <= candidate-limit <= 100; min-reactions >= 1; chebi-limit >= 0')
    weights = DEFAULT_WEIGHTS if weights is None else weights
    if not isinstance(weights, dict) or set(weights) != set(DEFAULT_WEIGHTS) or any(not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in weights.values()) or sum(weights.values()) <= 0:
        raise DiscoveryError('Weights must contain nonnegative coverage, continuity, evidence and central values')
    weights = {key: value / sum(weights.values()) for key, value in weights.items()}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    warnings, candidates, selected = [], [], []
    def write(name, data):
        (output / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    try:
        progress('Resolving organism from model annotations and referenced genes...', flush=True)
        taxon = infer_organism(model, client, organism, kegg_organism)
        write('organism.json', taxon)
        if taxon['status'] in ('unresolved', 'ambiguous'):
            raise DiscoveryError('Organism cannot be resolved; see organism.json and specify --organism')
        org = taxon.get('kegg_organism')
        warnings.extend(taxon['warnings'])
        mappings = compound_mapping(model)
        currency = automatic_currency(model, mappings)
        write('id_mapping.json', dict(metabolites=mappings, reactions=[dict(
            model_id=r['id'], kegg_reactions=annotation_ids(r, 'kegg.reaction'),
            ec_numbers=annotation_ids(r, 'ec-code'), gene_reaction_rule=r.get('gene_reaction_rule', ''))
            for r in model['reactions']], genes=[dict(model_id=g['id'], name=g.get('name', ''),
            ncbi_gene_ids=annotation_ids(g, 'ncbigene') + annotation_ids(g, 'ncbi.gene'),
            uniprot_ids=annotation_ids(g, 'uniprot'), locus_tags=annotation_ids(g, 'refseq_locus_tag'),
            kegg_genes=annotation_ids(g, 'kegg.genes')) for g in model.get('genes', [])]))
        if 'reactome' in databases:
            progress('Searching supplementary Reactome pathway membership...', flush=True)
            write('reactome_candidates.json', reactome_candidates(model, mappings, taxon, client, warnings, reactome_limit))
        progress('Searching KEGG pathway membership...', flush=True)
        pathway_names = {a.removeprefix('path:'): b for a, b in pairs(client.get(KEGG + 'list/pathway' + ('/' + org if org else '')))}
        compound_pathways, reaction_pathways = defaultdict(set), defaultdict(set)
        compounds = {m['kegg_compounds'][0] for m in mappings if m['status'] == 'unique' and m['model_id'] not in currency}
        for batch in batches(compounds):
            for a, b in pairs(client.get(KEGG + 'link/pathway/' + '+'.join('cpd:' + v for v in batch))):
                compound_pathways[b[-5:]].add(a.removeprefix('cpd:'))
        known_reactions = {kid for r in model['reactions'] for kid in annotation_ids(r, 'kegg.reaction') if re.fullmatch(r'R\d{5}', kid)}
        for batch in batches(known_reactions):
            for a, b in pairs(client.get(KEGG + 'link/pathway/' + '+'.join('rn:' + v for v in batch))):
                reaction_pathways[b[-5:]].add(a.removeprefix('rn:'))
        shortlist = []
        for pid, name in pathway_names.items():
            suffix = pid[-5:]
            # Overview maps are not individual representative pathways.
            if suffix.startswith(('011', '012', '013')) or suffix not in compound_pathways:
                continue
            count = len(compound_pathways[suffix])
            if count < 2:
                continue
            preliminary = count + 2 * len(reaction_pathways[suffix]) + (8 if suffix in CENTRAL else 0)
            shortlist.append((preliminary, pid, name))
        for _, pid, name in sorted(shortlist, key=lambda x: (-x[0], x[1]))[:candidate_limit]:
            progress('Matching KEGG pathway ' + pid + ': ' + name, flush=True)
            kgml_id = pid if org else 'ko' + pid[-5:]
            try:
                records = parse_kgml(client.get(KEGG + 'get/' + kgml_id + '/kgml'), org)
            except DiscoveryError as exc:
                warnings.append(str(exc)); continue
            matches = match_pathway(model, records, mappings, currency, org)
            db_ids = {r['kegg_reaction'] for r in records}
            # Only corroborated pairs enter the backbone; weak pairs remain in audit.
            corroborated = [m for m in matches if m['confidence'] >= 0.8]
            connection = continuity(corroborated)
            coverage = len({m['kegg_reaction'] for m in corroborated}) / max(len(db_ids), 1)
            evidence_score = sum(m['confidence'] for m in corroborated) / max(len(corroborated), 1)
            metrics = dict(coverage=coverage, continuity=connection['largest_component_fraction'],
                           evidence=evidence_score, central=float(pid[-5:] in CENTRAL))
            candidates.append(dict(id=pid, name=name, source='KEGG',
                                   source_url='https://www.kegg.jp/pathway/' + kgml_id,
                                   scope='organism_gene_supported' if org else 'reference_pathway',
                                   database_reaction_count=len(db_ids), matched_reaction_count=len(corroborated),
                                   metrics=metrics, score=sum(weights[k] * metrics[k] for k in weights),
                                   connectivity=connection, matches=matches, selected_steps=[]))
        candidates.sort(key=lambda c: (-c['score'], -c['matched_reaction_count'], c['id']))
        assigned = set()
        groups = []
        for candidate in candidates:
            eligible = [m for m in candidate['matches'] if m['confidence'] >= 0.8 and m['reaction_id'] not in assigned]
            parts = continuity(eligible)['reaction_components']
            biggest = set(sorted(parts, key=lambda p: (-len(p), p))[0]) if parts else set()
            steps = [m for m in eligible if m['reaction_id'] in biggest]
            if len(groups) < max_pathways and len(steps) >= min_reactions:
                candidate['selected_steps'] = steps
                assigned.update(s['reaction_id'] for s in steps)
                selected.append(candidate)
                groups.append(dict(id=candidate['id'], name=candidate['name'],
                                   steps=[dict(s) for s in steps], source_url=candidate['source_url']))
        cycles = backbone_layout(groups) if groups else []
        important = important_compounds(model, mappings, selected, cycles, currency)
        if 'chebi' in databases:
            progress('Retrieving ChEBI roles for important compound candidates...', flush=True)
            enrich_chebi(important, client, warnings, chebi_limit)
        write('important_metabolites.json', dict(metabolites=important,
              interpretation='Importance is contextual; roles of ambiguous ChEBI identities are candidate evidence'))
        write('pathway_candidates.json', dict(candidates=candidates, weights=weights, candidate_limit=candidate_limit,
              selection='Ranked candidates; largest connected component; reactions assigned once; minimum corroborated reaction count',
              limitations=['KGML pairs do not prove full stoichiometric identity.',
                           'Database coverage and structural connectivity do not establish flux or essentiality.',
                           'Ambiguous IDs require a unique KGML intersection plus independent gene/reaction evidence; global identities remain ambiguous.',
                           'Only shortlisted pathways are evaluated; ranking is not exhaustive.']))
        specification = dict(name=model['id'] + ': database-supported pathways', pathways=groups,
                             cycle_components=cycles, shared_primary_species=[],
                             selection_basis='KEGG reaction/gene evidence plus unique compound annotations and allowed model directions')
        write('auto_pathways.json', specification)
        if not groups:
            raise DiscoveryError('No sufficiently supported connected pathways; see pathway_candidates.json')
        return specification, currency
    finally:
        write('database_provenance.json', dict(requests=list(client.requests.values()), warnings=warnings,
                                             offline=client.offline, cache=str(client.cache)))
