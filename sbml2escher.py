#!/usr/bin/env python3
"""Discover bacterial central-carbon pathways and draw sampled or structural maps."""
import json
import copy
import re
import xml.etree.ElementTree as ET
import math
from pathlib import Path

import escher_grid as engine
from annotation_context import enrich_context
from bacterial_pathways import REFERENCE, discover, layout_plan
from model_input import load_model


def discovery_context(filename, model):
    """Read SBML compartment/SBO context without modifying the sampled model hash."""
    if Path(filename).suffix.lower()=='.json':
        return model
    context=copy.deepcopy(model);root=ET.parse(filename).getroot()
    def local(tag):
        return tag.rsplit('}',1)[-1]
    def decode(value):
        return re.sub(r'__(\d+)__',lambda match:chr(int(match[1])),value.removeprefix('R_'))
    context['compartments']={e.get('id'):e.get('name',e.get('id')) for e in root.iter() if local(e.tag)=='compartment'}
    reactions={r['id']:r for r in context['reactions']}
    for entity in root.iter():
        if local(entity.tag)=='reaction' and entity.get('sboTerm') and decode(entity.get('id','')) in reactions:
            reactions[decode(entity.get('id'))]['sbo_term']=entity.get('sboTerm')
    context['_unsupported_sampling_constraints'] = any(local(e.tag) == 'userDefinedConstraint' for e in root.iter())
    context['_boundary_species']=[e.get('id') for e in root.iter() if local(e.tag)=='species' and e.get('boundaryCondition') in ('true','1')]
    return enrich_context(context, root)


class BacterialStrategy:
    def __init__(self, model, reference=REFERENCE, identity_map=None, glucose_id=None, max_milp_nodes=180):
        self.context = model
        self.routes,self.extraction,self.mapped,self.reference=discover(model,reference,identity_map,glucose_id)
        self.glucose_id=self.extraction['glucose_metabolite_id']
        self.max_milp_nodes=max_milp_nodes
        for route in self.routes:
            route['model_roles']={mid:self.mapped[mid] for mid in route['metabolite_ids']}

    def classify(self, model, overrides=None):
        classes=engine.classify_metabolites(self.context)
        for mid,role in self.mapped.items():
            classes[mid]['kind']=self.reference['roles'][role]['kind']
            classes[mid]['reason']='Mapped BiGG/KEGG chemical role: '+role
        for metabolite in self.context['metabolites']:
            mid = metabolite['id']
            if metabolite.get('_boundary_placeholder'):
                classes[mid].update(boundary_placeholder=True, reason=metabolite['_boundary_evidence'])
            if metabolite.get('_formula_evidence'):
                classes[mid]['formula_evidence'] = metabolite['_formula_evidence']
        if overrides:
            engine.classify_metabolites(model,overrides)  # Validate IDs and allowed values.
            for mid,kind in overrides.items():
                classes[mid].update(kind=kind,reason='explicit classification override')
        return classes

    def pathways(self, model):
        return self.routes,self.extraction['missing_pathways']

    def exchanges(self, grid_ids):
        return {'rxn:'+e['reaction_id'] for e in self.extraction['exchanges'] if 'rxn:'+e['reaction_id'] in grid_ids}

    def plan(self, routes, grid_ids):
        # Column zero is now reserved for exchanges. Keep the entire central
        # skeleton collinear/square one column inside the perimeter.
        fixed, plan = layout_plan(routes,grid_ids)
        shift_y = -int(not plan['glucose_top_margin'])
        fixed = {key:[x+1,y+shift_y] for key,(x,y) in fixed.items()}
        plan['fixed_node_offsets'] = fixed
        for route in plan['routes']:
            route['node_offsets'] = {key:[x+1,y+shift_y] for key,(x,y) in route['node_offsets'].items()}
        if plan['square']:
            plan['square']['left'] += 1
            plan['square']['top_offset'] += shift_y
        plan['glucose_top_margin'] = max(1,plan['glucose_top_margin'])
        return fixed, plan

    def grid_options(self, nodes, backbones):
        margin=backbones['glucose_top_margin']
        exchange_ids={e['reaction_id'] for e in self.extraction['exchanges']}
        return dict(top_left_margin=max(1,margin), exclusive_boundary=True,
                    node_top_margins={key:margin for key,n in nodes.items() if n['kind']!='currency' and
                                      (n['kind']!='reaction' or n['original_id'] not in exchange_ids)},
                    max_milp_nodes=self.max_milp_nodes)

    def save(self, output):
        engine.dump(output/'pathway_extraction.json',self.extraction)
        engine.dump(output/'id_mapping.json',self.extraction['mapping'])


def parser():
    p=engine.parser()
    p.description=__doc__
    p.prog='sbml2escher.py'
    p.set_defaults(output_dir=Path(__file__).parent/'e_coli_core_sbml2escher_output',program_name='sbml2escher')
    p.add_argument('--reference',type=Path,default=REFERENCE,help='Versioned BiGG/KEGG compound-role reference JSON')
    p.add_argument('--identity-map',type=Path,help='JSON mapping model-local metabolite IDs to reference chemical roles')
    p.add_argument('--discover-only',action='store_true',help='Extract and report pathways without sampling or layout')
    p.add_argument('--max-milp-nodes',type=int,default=180,help='Above this grid-node count, preserve constraints with heuristic placement')
    return p


def run(args):
    if not math.isfinite(args.coverage) or not 0<=args.coverage<=100:
        raise ValueError('coverage must be in [0,100]')
    if not math.isfinite(args.time_limit) or args.time_limit<=0 or args.max_milp_nodes<1:
        raise ValueError('time-limit and max-milp-nodes must be positive')
    if not args.no_sampling and (args.samples<20 or args.chains<2 or args.thinning<1 or args.burn_in<1):
        raise ValueError('Require samples >=20, chains >=2, thinning >=1, burn-in >=1')
    if not math.isfinite(args.zero_tolerance) or args.zero_tolerance<=0:
        raise ValueError('zero-tolerance must be positive and finite')
    if any(v is not None and v<1 for v in (args.grid_cols,args.grid_rows)):
        raise ValueError('Grid dimensions must be positive')
    engine.validate_sampling_options(args)
    model=load_model(args.model)
    identity=json.loads(args.identity_map.read_text()) if args.identity_map else None
    context=discovery_context(args.model,model)
    strategy=BacterialStrategy(context,args.reference,identity,args.glucose_id,args.max_milp_nodes)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if args.discover_only:
        strategy.save(args.output_dir)
        print(json.dumps(dict(output=str(args.output_dir),status=strategy.extraction['status'],
                              pathways=[dict(id=r['id'],reactions=len(r['reaction_ids']),closed=r['closed']) for r in strategy.routes]),indent=2))
        return strategy.extraction
    if context.get('_boundary_species') and not args.no_sampling:
        args._sampling_unavailable_reason = 'Sampling boundaryCondition=true SBML species is unsupported'
    if context.get('_unsupported_sampling_constraints') and not args.no_sampling:
        args._sampling_unavailable_reason = 'SBML userDefinedConstraint is not supported by flux sampling'
    if any(m.get('_boundary_placeholder') for m in context['metabolites']) and not args.no_sampling:
        args._sampling_unavailable_reason = 'CellDesigner source/sink placeholders require explicit mass-balance semantics; structural display only'
    args._strategy=strategy
    return engine.run(args)


def main():
    p=parser();args=p.parse_args()
    try:
        run(args)
    except (ValueError,RuntimeError,OSError) as exc:
        p.error(str(exc))


if __name__=='__main__':
    main()
