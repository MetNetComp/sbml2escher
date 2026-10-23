"""Semantic fallback tests using small SBML fixtures, with no network access."""
import copy
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from annotation_context import CD, hidden_metabolite
from bacterial_pathways import REFERENCE, discover, map_roles
from escher_grid import display_graph, select_structural_reactions
from model_input import load_model
from sbml2escher import BacterialStrategy, discovery_context, parser, run

CORE = '{http://www.sbml.org/sbml/level3/version1/core}'


def fixture():
    root = ET.Element(CORE+'sbml', level='3', version='1')
    model = ET.SubElement(root, CORE+'model', id='opaque_model')
    compartments = ET.SubElement(model, CORE+'listOfCompartments')
    for cid, name in [('default', 'Default'), ('c2', 'Extracellular space'), ('c3', 'Cytoplasm')]:
        ET.SubElement(compartments, CORE+'compartment', id=cid, name=name, constant='true')
    species = ET.SubElement(model, CORE+'listOfSpecies')
    labels = ['glc_D[e]', 'glc_D[c]', 'g6p[c]', '6pgl[c]', '6pgc[c]', '2ddg6p[c]', 'pyr[c]', 'accoa[c]', 'atp[c]', 'adp[c]', 'h[c]', '.']
    for i, label in enumerate(labels):
        e = ET.SubElement(species, CORE+'species', id=f'x{i}', name='.' if label=='.' else f'opaque_{i}',
                          compartment='default' if label=='.' else 'c2' if i==0 else 'c3')
        annotation = ET.SubElement(e, CORE+'annotation'); extension = ET.SubElement(annotation, CD+'extension')
        identity = ET.SubElement(extension, CD+'speciesIdentity')
        ET.SubElement(identity, CD+'class').text = 'UNKNOWN' if label=='.' else 'SIMPLE_MOLECULE'
        ET.SubElement(identity, CD+'name').text = label
        ET.SubElement(extension, CD+'bounds', x='999999', y='-999999')
    reactions = ET.SubElement(model, CORE+'listOfReactions')
    def reaction(rid, label, coefficients):
        e = ET.SubElement(reactions, CORE+'reaction', id=rid, name='', reversible='true')
        extension = ET.SubElement(ET.SubElement(e, CORE+'annotation'), CD+'extension')
        ET.SubElement(extension, CD+'name').text = label
        ET.SubElement(extension, CD+'reactionType').text = 'TRANSPORT'
        ET.SubElement(extension, CD+'line', color='ffff0000')
        for group, tag, sign in [('listOfReactants','baseReactant',-1),('listOfProducts','baseProduct',1)]:
            refs = ET.SubElement(e, CORE+group)
            for mid, coefficient in coefficients.items():
                if coefficient*sign>0:
                    ET.SubElement(refs,CORE+'speciesReference',species=mid,stoichiometry=str(abs(coefficient)),constant='true')
                    ET.SubElement(extension,CD+tag,species=mid)
        return e
    reaction('exchange', 'EX_glc_D[e]', {'x0':-1,'x11':1})
    for i in range(7):
        coefficients={f'x{i}':-1,f'x{i+1}':1}
        if i==1:coefficients.update(x8=-1,x9=1,x10=1)
        reaction(f'enzyme{i}',f'enzyme_{i}',coefficients)
    return root


class AnnotationContextTests(unittest.TestCase):
    def context(self, root):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        path=Path(temporary.name)/'fixture.xml';ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
        source=load_model(path);before=copy.deepcopy(source);context=discovery_context(path,source)
        self.assertEqual(source,before)
        return source,context

    def test_cd_only_labels_recover_real_directed_path_without_source_changes(self):
        source,context=self.context(fixture());strategy=BacterialStrategy(context)
        self.assertEqual([r['id'] for r in strategy.routes],['glycolysis','pdh','glucose_exchange'])
        primary=strategy.routes[0];self.assertEqual(primary['family'],'ED')
        self.assertEqual(primary['reaction_ids'],[f'enzyme{i}' for i in range(6)])
        self.assertEqual(primary['exchange_reaction_id'],'exchange')
        record=next(x for x in strategy.extraction['mapping'] if x['model_id']=='x0')
        self.assertEqual(record['role'],'glc__D')
        self.assertTrue(any(e['source']=='CellDesigner speciesIdentity/name' for e in record['evidence']))
        self.assertFalse(strategy.extraction['annotation_fallback']['drawing_coordinates_used'])
        self.assertFalse(any('_identity_hints' in m for m in source['metabolites']))
        self.assertEqual([(r['metabolites'],r['lower_bound'],r['upper_bound']) for r in source['reactions']],
                         [(r['metabolites'],r['lower_bound'],r['upper_bound']) for r in context['reactions']])

    def test_boundaries_hidden_only_in_projection_and_cofactors_classified(self):
        source,context=self.context(fixture());strategy=BacterialStrategy(context);classes=strategy.classify(source)
        self.assertTrue(classes['x11']['boundary_placeholder']);self.assertEqual(classes['x11']['kind'],'organic')
        self.assertEqual(classes['x8']['kind'],'coenzyme');self.assertEqual(classes['x10']['kind'],'inorganic')
        selected,records,_=select_structural_reactions(source,classes,100)
        projection,nodes,edges=display_graph(source,selected,classes)
        self.assertNotIn('met:x11',nodes);self.assertNotIn('met:x10',nodes)
        exchange=next(x for x in projection['reactions'] if x['id']=='exchange')
        self.assertEqual(exchange['metabolites'],{'x0':-1})
        self.assertEqual(source['reactions'][0]['metabolites'],{'x0':-1,'x11':1})
        record=next(x for x in records if x['reaction_id']=='exchange')
        self.assertEqual(record['omitted_boundary_placeholders'],['x11']);self.assertEqual(record['omitted_inorganic'],[])

    def test_structured_notes_and_vmh_link_are_used_without_cd_names(self):
        root=fixture();species=root.find('.//'+CORE+'species')
        identity=species.find('.//'+CD+'speciesIdentity');identity.remove(identity.find(CD+'name'))
        notes=ET.SubElement(species,CORE+'notes');notes.text='Name: D-glucose Abbreviation: glc_D ChargedFormula: C6H12O6 Charge: 0 '
        ET.SubElement(notes,'{http://www.w3.org/1999/xhtml}a',href='https://delta.vmh.life/metabolite/glc_D').text='VMH Metabolite'
        source,context=self.context(root);mapped,records=map_roles(context,json.loads(REFERENCE.read_text()))
        self.assertEqual(mapped['x0'],'glc__D');self.assertEqual(context['metabolites'][0]['chemicalFormula'],'C6H12O6')
        self.assertNotIn('chemicalFormula',source['metabolites'][0])
        evidence=records[0]['evidence'];self.assertTrue(any(x['source']=='SBML notes:VMH link' and x.get('uri') for x in evidence))

    def test_conflicting_identifiers_stay_ambiguous_and_override_resolves(self):
        root=fixture();species=root.find('.//'+CORE+'species');ET.SubElement(species,CORE+'notes').text='Abbreviation: pyr'
        source,context=self.context(root);reference=json.loads(REFERENCE.read_text());mapped,records=map_roles(context,reference)
        self.assertNotIn('x0',mapped);self.assertEqual(records[0]['status'],'ambiguous')
        self.assertEqual(set(records[0]['candidate_roles']),{'pyr','glc__D'})
        mapped,records=map_roles(context,reference,{'x0':'glc__D'});self.assertEqual(mapped['x0'],'glc__D')
        context['metabolites'][0]['annotation']={'kegg.compound':['C00022']}
        mapped,_=map_roles(context,reference);self.assertEqual(mapped['x0'],'pyr')

    def test_transport_type_or_paint_alone_never_identifies_exchange_or_pathway(self):
        root=fixture()
        for e in root.iter(CD+'name'):
            e.text='unknown_label'
        source,context=self.context(root);routes,report,mapped,_=discover(context)
        self.assertFalse(routes);self.assertFalse(report['exchanges']);self.assertFalse(mapped)
        self.assertFalse(context['metabolites'][-1].get('_boundary_placeholder'))

    def test_wrong_cd_connections_do_not_establish_exchange(self):
        root=fixture();reaction=root.find('.//'+CORE+'reaction')
        reaction.find('.//'+CD+'baseReactant').set('species','not_in_source_reaction')
        source,context=self.context(root);_,report,_,_=discover(context)
        self.assertFalse(report['exchanges'])

    def test_complete_cli_pipeline_preserves_source_and_falls_back_for_symbols(self):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        base=Path(temporary.name);path=base/'fixture.xml';out=base/'output'
        ET.ElementTree(fixture()).write(path,encoding='utf-8',xml_declaration=True)
        args=parser().parse_args(['--model',str(path),'--output-dir',str(out),'--max-milp-nodes','1'])
        with patch('escher_grid.sample_fluxes',side_effect=AssertionError('source/sink symbols must not be sampled')), contextlib.redirect_stdout(io.StringIO()):
            report=run(args)
        self.assertEqual(report['sampling']['mode'],'fallback')
        self.assertIn('placeholders',report['sampling']['reason'])
        self.assertTrue(report['validation']['escher']['schema_valid'])
        self.assertTrue(report['validation']['no_boundary_placeholder_nodes'])
        self.assertGreater(report['selection']['major_selected_reactions'],0)
        self.assertTrue(report['exchange_connectors']['all_exchanges_connected'])
        self.assertEqual(load_model(path),json.loads((out/'source_model.json').read_text()))
        self.assertFalse((out/'flux_statistics.json').exists())
        self.assertTrue((out/'index.html').is_file())

    def test_protein_label_is_not_a_chemical_identity(self):
        root=fixture();root.find('.//'+CD+'class').text='PROTEIN'
        _,context=self.context(root);mapped,_=map_roles(context,json.loads(REFERENCE.read_text()))
        self.assertNotIn('x0',mapped)


if __name__=='__main__':unittest.main()
