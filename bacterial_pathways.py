"""Model-local, evidence-audited bacterial central-carbon backbone discovery.

BiGG/KEGG/ChEBI/MetaNetX/SEED annotations identify chemical roles, never replace
source reactions. Directed model incidences must realize every template step.
Missing or ambiguous evidence stays explicit; there is no universal guarantee.
"""
from collections import defaultdict
import hashlib
import heapq
import itertools
import json
import math
from pathlib import Path
import re

from auto_pathways import annotation_ids
from annotation_context import fallback_role_evidence, reaction_labels
from major_backbones2 import RED, validate_backbones

REFERENCE = Path(__file__).parent / 'data/bacteria/central_reference.json'
NAMESPACES = ('bigg.metabolite', 'kegg.compound', 'chebi', 'metanetx.chemical', 'seed.compound', 'biocyc', 'inchi_key')
EMP = ['glc__D', 'g6p', 'f6p', 'fdp', 'dhap', 'g3p', '13dpg', '3pg', '2pg', 'pep', 'pyr']
TCA = ['oaa', 'cit', 'acon_C', 'icit', 'akg', 'succoa', 'succ', 'fum', 'mal__L', 'oaa']


def normalized_name(value):
    value = re.sub(r'\[[^]]+\]$', '', value.strip())
    return re.sub(r'[^a-z0-9]', '', value.lower())


def map_roles(model, reference, overrides=None):
    roles = reference['roles']; index = defaultdict(set); names = defaultdict(set)
    for role, record in roles.items():
        for ns, values in record['aliases'].items():
            for value in values:
                index[(ns, value)].add(role)
        for name in record['names']:
            names[normalized_name(name)].add(role)
    overrides = overrides or {}
    known_ids = {m['id'] for m in model['metabolites']}
    if set(overrides) - known_ids or set(overrides.values()) - set(roles):
        raise ValueError('Identity map must assign known source metabolite IDs to roles in the reference catalog')
    mapped, records = {}, []
    for metabolite in model['metabolites']:
        mid = metabolite['id']; evidence = []
        comp = metabolite.get('compartment', '')
        base = mid.removesuffix('_' + comp) if comp else re.sub(r'_[a-z][0-9]?$', '', mid)
        namespace='__'+str(model.get('id',''))+'__'
        if base.startswith(namespace):
            base=base[len(namespace):]
        # BiGG SBML ID normalization is performed by model_input, not guessed here.
        if base in roles:
            evidence.append(dict(role=base, tier=300, source='BiGG-compatible metabolite ID', value=base))
        if re.fullmatch(r'cpd\d+', base):
            for role in index.get(('seed.compound', base), []):
                evidence.append(dict(role=role, tier=250, source='ModelSEED compound ID', value=base))
        for ns in NAMESPACES:
            values = annotation_ids(metabolite, ns)
            for value in values:
                for role in index.get((ns, value), []):
                    evidence.append(dict(role=role, tier=300 if ns=='bigg.metabolite' else 200 if ns in ('kegg.compound','metanetx.chemical','seed.compound','inchi_key') else 150,
                                         source='annotation:' + ns, value=value))
        # Exact reference names are a last resort; never fuzzy substring matching.
        for role in names.get(normalized_name(metabolite.get('name', '')), []):
            evidence.append(dict(role=role, tier=50, source='exact reference name', value=metabolite.get('name')))
        if not evidence:
            evidence.extend(fallback_role_evidence(metabolite, reference))
        if mid in overrides:
            evidence.append(dict(role=overrides[mid], tier=1000, source='explicit identity map', value=overrides[mid]))
        best = max((e['tier'] for e in evidence), default=0)
        candidates = sorted({e['role'] for e in evidence if e['tier']==best})
        status = 'mapped' if len(candidates)==1 else 'ambiguous' if candidates else 'unmapped'
        if status=='mapped':
            mapped[mid] = candidates[0]
        records.append(dict(model_id=mid, compartment=comp, role=mapped.get(mid), status=status,
                            candidate_roles=candidates, evidence=evidence,
                            lower_tier_conflicts=sorted({e['role'] for e in evidence} - set(candidates))))
        if metabolite.get('_identity_hints'):
            records[-1]['fallback_identity_hints'] = metabolite['_identity_hints']
            records[-1]['mapping_reason'] = ('recognized boundary placeholder, not a chemical role' if metabolite.get('_boundary_placeholder') else
                                            'ambiguous best evidence; explicit identity map required' if status == 'ambiguous' else
                                            'no exact reference match' if status == 'unmapped' else
                                            'highest-priority unambiguous evidence')
    return mapped, records


def import_allowed(reaction, mid):
    coefficient = reaction['metabolites'][mid]
    return reaction['lower_bound'] < 0 if coefficient < 0 else reaction['upper_bound'] > 0


def exchange_reactions(model, mapped):
    compounds = {m['id']: m for m in model['metabolites']}
    result = []
    for reaction in model['reactions']:
        if len(reaction['metabolites']) != 1:
            placeholders = [m for m in reaction['metabolites'] if compounds[m].get('_boundary_placeholder')]
            physical = [m for m in reaction['metabolites'] if m not in placeholders]
            labels = reaction_labels(reaction)
            explicit = [(source, value) for source, value in labels if value.startswith('EX_')]
            conflicting = any(value.startswith(('DM_', 'SK_')) for source, value in labels)
            if (len(placeholders) == 1 and len(physical) == 1 and explicit and not conflicting
                    and reaction.get('_celldesigner_connections_verified', True)
                    and reaction['metabolites'][placeholders[0]] * reaction['metabolites'][physical[0]] < 0):
                mid = physical[0]
                result.append(dict(reaction_id=reaction['id'], metabolite_id=mid,
                                   import_allowed=import_allowed(reaction, mid), boundary_kind='exchange',
                                   boundary_placeholder_id=placeholders[0],
                                   evidence='Explicit EX name plus verified CellDesigner boundary placeholder',
                                   name_evidence=[dict(source=s, value=v) for s, v in explicit]))
            continue
        mid = next(iter(reaction['metabolites'])); metabolite = compounds[mid]
        sbo = set(annotation_ids(reaction, 'sbo')) | {str(reaction.get('sbo_term', ''))}
        compartment=metabolite.get('compartment','')
        compartment_name=model.get('compartments',{}).get(compartment,compartment)
        external = compartment.lower() in ('e','e0','extracellular','external','outside','extracellular_space') or any(word in str(compartment_name).lower() for word in ('extracellular','external','outside'))
        explicit = reaction['id'].startswith('EX_') or bool({'SBO:0000627','0000627'} & sbo)
        if explicit or external or (mapped.get(mid)=='glc__D' and import_allowed(reaction, mid)):
            result.append(dict(reaction_id=reaction['id'], metabolite_id=mid,
                               import_allowed=import_allowed(reaction, mid),
                               evidence='EX identifier/SBO' if explicit else 'single external species' if external else 'single glucose source',
                               boundary_kind='exchange_candidate' if not explicit else 'exchange'))
    return result


def directed_graph(model, mapped, permitted_kinds=None, reference=None):
    graph = defaultdict(list)
    compounds = {m['id']: m for m in model['metabolites']}
    allowed = set(mapped)
    if permitted_kinds:
        allowed = {m for m, role in mapped.items() if reference['roles'][role]['kind'] in permitted_kinds}
    for reaction in model['reactions']:
        rid = reaction['id']
        if len(reaction['metabolites']) > 30 or 'biomass' in (rid + ' ' + reaction.get('name','')).lower():
            continue
        organic_roles = set()
        for mid in reaction['metabolites']:
            role=mapped.get(mid)
            if compounds[mid].get('_boundary_placeholder'):
                continue
            if role:
                if reference['roles'][role]['kind']=='organic':
                    organic_roles.add(role)
            else:
                formula=compounds[mid].get('chemicalFormula',compounds[mid].get('formula','')) or ''
                if not formula or 'C' in re.findall(r'[A-Z][a-z]?',formula):
                    organic_roles.add('unmapped:'+mid)
        substrates = [m for m,c in reaction['metabolites'].items() if c<0 and m in allowed]
        products = [m for m,c in reaction['metabolites'].items() if c>0 and m in allowed]
        for left, right, direction in [(substrates,products,1),(products,substrates,-1)]:
            if (direction==1 and reaction['upper_bound']<=0) or (direction==-1 and reaction['lower_bound']>=0):
                continue
            for a in left:
                for b in right:
                    preferred=reference.get('reaction_preferences',{}).get(mapped[a]+'|'+mapped[b],[])
                    reaction_ids=annotation_ids(reaction,'kegg.reaction')
                    preference=0 if not preferred or set(preferred)&set(reaction_ids) else .01
                    graph[a].append(dict(source=a,target=b,reaction_id=rid,direction=direction,
                                         organic_roles=sorted(organic_roles),reference_penalty=preference))
    for edges in graph.values():
        edges.sort(key=lambda e:(e['reaction_id'],e['target']))
    return graph


def supported_carbon_context(edge, mapped):
    """Exclude PTS PEP->pyruvate side pairs and unrelated CoA transfer reactions."""
    a,b=mapped[edge['source']],mapped[edge['target']]
    allowed={a,b}
    if (a,b)==('glc__D','g6p'):
        allowed.update(['pep','pyr'])
    if (a,b) in (('fdp','dhap'),('fdp','g3p')):
        allowed.update(['dhap','g3p'])
    if (a,b) in (('2ddg6p','pyr'),('2ddglcn','pyr')):
        allowed.update(['g3p','glyc'])
    if (a,b)==('pyr','accoa'):
        allowed.add('for')
    if (a,b)==('oaa','cit'):
        allowed.add('accoa')
    return set(edge.get('organic_roles',[])) <= allowed


def match_pattern(graph, mapped, pattern, starts=None, closed=False, max_states=30000, diagnostics=None):
    """Search ordered chemical roles; same-role transport/anomer steps are explicit."""
    starts = sorted(starts if starts is not None else [m for m,r in mapped.items() if r==pattern[0]])
    heap = []; serial = itertools.count()
    for start in starts:
        if mapped.get(start)==pattern[0]:
            heapq.heappush(heap,(0,next(serial),start,0,[start],[]))
    states = 0; seen = {}
    while heap:
        cost, _, mid, stage, species, steps = heapq.heappop(heap); states+=1
        if states>max_states:
            if diagnostics is not None:
                diagnostics.append(dict(pattern=pattern,status='search_limit_reached',max_states=max_states))
            return None
        if stage==len(pattern)-1:
            if not closed or species[-1]==species[0]:
                return dict(metabolite_ids=species,reaction_ids=[s['reaction_id'] for s in steps],steps=steps,
                            pattern=pattern,search_states=states)
            continue
        used = frozenset(s['reaction_id'] for s in steps)
        signature=(mid,stage,used)
        if seen.get(signature,math.inf)<cost:
            continue
        seen[signature]=cost
        for edge in graph.get(mid,[]):
            if not supported_carbon_context(edge,mapped):
                continue
            target=edge['target']; role=mapped.get(target)
            if edge['reaction_id'] in used:
                continue
            if role==pattern[stage+1]:
                nxt=stage+1
            elif role==pattern[stage]:
                nxt=stage
            else:
                continue
            returning=closed and nxt==len(pattern)-1 and target==species[0]
            if target in species and not returning:
                continue
            heapq.heappush(heap,(cost+1+edge.get('reference_penalty',0),next(serial),target,nxt,species+[target],steps+[edge]))
    return None


def discover(model, reference_path=REFERENCE, identity_map=None, glucose_id=None):
    reference_path=Path(reference_path); raw=reference_path.read_bytes(); reference=json.loads(raw)
    overrides=dict(identity_map or {})
    if glucose_id:
        overrides[glucose_id]='glc__D'
    mapped, mapping=map_roles(model,reference,overrides)
    graph=directed_graph(model,mapped,{'organic'},reference)
    exchanges=exchange_reactions(model,mapped)
    glucose_exchanges=[e for e in exchanges if mapped.get(e['metabolite_id'])=='glc__D']
    glucose_starts=[e['metabolite_id'] for e in glucose_exchanges if e['import_allowed']]
    if not glucose_starts:
        glucose_starts=[e['metabolite_id'] for e in glucose_exchanges]
    if not glucose_starts:
        glucose_starts=[m for m,r in mapped.items() if r=='glc__D']
    candidates=[];search_warnings=[]
    patterns=[('EMP',EMP),('EMP_direct_G3P',[x for x in EMP if x!='dhap']),
              ('EMP_GAPN',[x for x in EMP if x not in ('dhap','13dpg')]),
              ('ED',['glc__D','g6p','6pgl','6pgc','2ddg6p','pyr']),
              ('ED_gluconate',['glc__D','glcn','6pgc','2ddg6p','pyr']),
              ('ED_semiphosphorylative',['glc__D','glcn','2ddglcn','2ddg6p','pyr']),
              ('ED_nonphosphorylative',['glc__D','glcn','2ddglcn','pyr'])]
    for family, pattern in patterns:
        match=match_pattern(graph,mapped,pattern,glucose_starts,diagnostics=search_warnings)
        if match:
            candidates.append(dict(family=family,**match))
    routes=[]; missing=[]; chosen_glucose=None
    if candidates:
        primary=candidates[0]
        chosen_glucose=primary['metabolite_ids'][0]
        route=dict(primary,id='glycolysis',name='Glucose to pyruvate ('+primary['family']+')',shape='vertical',closed=False)
        options=sorted([e for e in glucose_exchanges if e['metabolite_id']==chosen_glucose],key=lambda e:(not e['import_allowed'],e['reaction_id']))
        route['exchange_reaction_id']=options[0]['reaction_id'] if options else None
        routes.append(route)
        start_pyr=primary['metabolite_ids'][-1]
    else:
        missing.append(dict(pathway='glucose_to_pyruvate',reason='No complete directed EMP/ED pattern in the supplied model',
                            identified_glucose_species=glucose_starts,identified_pyruvate_species=[m for m,r in mapped.items() if r=='pyr']))
        chosen_glucose=glucose_starts[0] if glucose_starts else None
        start_pyr=None
    # Pyruvate oxidation alternatives are selected by source connectivity, not PDH IDs.
    pyr_match=None
    for pattern in [['pyr','accoa'],['pyr','hethmpp','acdhlam','accoa']]:
        pyr_match=match_pattern(graph,mapped,pattern,[start_pyr] if start_pyr else None,diagnostics=search_warnings)
        if pyr_match:
            break
    if pyr_match:
        routes.append(dict(pyr_match,id='pdh',name='Pyruvate to acetyl-CoA (model reaction)',shape='horizontal',closed=False))
    else:
        missing.append(dict(pathway='pyruvate_to_acetyl_coa',reason='No directed model-local pyruvate/acetyl-CoA connection'))
    # An oxidative TCA cycle must return to the same model-local oxaloacetate species.
    cycle=None
    tca_multistep=TCA[:5]+['sucdhlam']+TCA[5:]
    for pattern in [TCA,[x for x in TCA if x!='acon_C'],tca_multistep,[x for x in tca_multistep if x!='acon_C']]:
        cycle=match_pattern(graph,mapped,pattern,closed=True,diagnostics=search_warnings)
        if cycle:
            break
    if cycle:
        routes.append(dict(cycle,id='tca',name='Directed model-local TCA cycle',shape='square',closed=True))
    else:
        partial=[]
        for a,b in zip(TCA,TCA[1:]):
            segment=match_pattern(graph,mapped,[a,b],diagnostics=search_warnings)
            if segment:
                partial.extend(segment['steps'])
        if not any(mapped[s['source']]=='cit' and mapped[s['target']]=='acon_C' for s in partial):
            direct=match_pattern(graph,mapped,['cit','icit'],diagnostics=search_warnings)
            if direct:
                partial.extend(direct['steps'])
        if partial:
            # Do not manufacture a closure across absent model reactions.
            used_mids=list(dict.fromkeys([s['source'] for s in partial]+[s['target'] for s in partial]))
            routes.append(dict(id='tca',name='Partial TCA components (not a complete cycle)',shape='square',closed=False,
                               metabolite_ids=used_mids,reaction_ids=list(dict.fromkeys(s['reaction_id'] for s in partial)),
                               steps=partial,pattern=TCA))
        missing.append(dict(pathway='tca',reason='No complete directed oxidative cycle; available segments retained without invented edges',
                            available_steps=partial))
    if chosen_glucose:
        options=sorted([e for e in glucose_exchanges if e['metabolite_id']==chosen_glucose],key=lambda e:(not e['import_allowed'],e['reaction_id']))
        if options:
            routes.append(dict(id='glucose_exchange',name='Glucose exchange',shape='exchange',metabolite_ids=[chosen_glucose],
                               reaction_ids=[options[0]['reaction_id']],steps=[],closed=False))
    for route in routes:
        route['evidence']='Reference chemical-role pattern + actual directed source-model incidences'
        route['flux_feasibility_proven']=False
    extraction=dict(model_id=model.get('id'),reference_sha256=hashlib.sha256(raw).hexdigest(),reference_sources=reference['sources'],
                    reference_modules=reference['modules'],mapping=mapping,mapped_metabolites=len(mapped),
                    total_metabolites=len(model['metabolites']),glucose_metabolite_id=chosen_glucose,
                    exchanges=exchanges,glucose_candidates=candidates,routes=routes,missing_pathways=missing,search_warnings=search_warnings,
                    status='complete_central_backbone' if len([r for r in routes if r['id'] in ('glycolysis','pdh','tca') and (r['id']!='tca' or r['closed'])])==3 else 'partial_or_unresolved',
                    limitations=['Chemical role mapping does not assert identical stereochemistry/protonation.',
                                 'Connectivity and bounds do not prove simultaneous physiological activity or carbon atom transfer.',
                                 'Absent or ambiguous annotations and genuinely missing pathways cannot be recovered with certainty.',
                                 'EMP is preferred when complete; otherwise ED variants are tried in documented order.'])
    fallback_records = [r for r in mapping if any(e['source'].startswith(('SBML notes:', 'SBML species name', 'CellDesigner ')) for e in r['evidence'])]
    placeholders = [m['id'] for m in model['metabolites'] if m.get('_boundary_placeholder')]
    if fallback_records or placeholders:
        extraction['annotation_fallback'] = dict(
            mapped_with_supplemental_evidence=sum(r['status'] == 'mapped' for r in fallback_records),
            ambiguous_with_supplemental_evidence=sum(r['status'] == 'ambiguous' for r in fallback_records),
            boundary_placeholder_ids=placeholders,
            exchange_reactions_with_boundary_placeholders=[e['reaction_id'] for e in exchanges if e.get('boundary_placeholder_id')],
            source_chemistry_modified=False, drawing_coordinates_used=False)
    return routes,extraction,mapped,reference


def layout_plan(routes, displayed_nodes):
    """EX -> glucose -> ... -> pyruvate is one ordered vertical grid line."""
    offsets={}; per_route={}; square=None
    def put(key,point,local):
        point=list(point)
        if key in offsets and offsets[key]!=point:
            raise ValueError('Overlapping pathway constraints disagree: '+key)
        offsets[key]=point;local[key]=point
    ordered=sorted(routes,key=lambda r: {'glycolysis':0,'glucose_exchange':1,'pdh':2,'tca':3}.get(r['id'],4))
    last_acetyl=None; margin=0
    for route in ordered:
        if route['id']=='tca':
            continue
        mids=route['metabolite_ids']; local={}
        if route['id']=='glycolysis':
            exchange=route.get('exchange_reaction_id');margin=int(bool(exchange))
            chain=(['rxn:'+exchange] if exchange else [])
            for mid,rid in zip(mids,route['reaction_ids']):
                chain.extend(['met:'+mid,'rxn:'+rid])
            chain.append('met:'+mids[-1])
            for i,key in enumerate(chain):
                put(key,[0,-i],local)
        elif route['id']=='glucose_exchange':
            margin=1
            put('rxn:'+route['reaction_ids'][0],[0,0],local);put('met:'+mids[0],[0,-1],local)
        else:
            start=offsets.get('met:'+mids[0],[0,-2-margin])
            chain=[]
            for mid,rid in zip(mids,route['reaction_ids']):
                chain.extend(['met:'+mid,'rxn:'+rid])
            chain.append('met:'+mids[-1])
            for i,key in enumerate(chain):
                put(key,[start[0]+i,start[1]],local)
            last_acetyl=offsets['met:'+mids[-1]]
        per_route[route['id']]=local
    for route in ordered:
        if route['id']!='tca':
            continue
        mids=route['metabolite_ids'];steps=route['steps'];local={}
        if route['closed']:
            mids=mids[:-1]
            chain=[key for mid,step in zip(mids,steps) for key in ('met:'+mid,'rxn:'+step['reaction_id'])]
        else:
            # Components keep their supported species/reactions. Missing steps leave gaps.
            roles=route['pattern'][:-1];by_role={}
            for mid in mids:
                # Filled by strategy with an explicit role map, never inferred from IDs.
                role=route['model_roles'][mid];by_role.setdefault(role,mid)
            chain=[]
            for i,role in enumerate(roles):
                mid=by_role.get(role)
                targets={route['pattern'][i+1]}
                if role=='cit':
                    targets.add('icit')
                step=next((s for s in steps if s['source']==mid and route['model_roles'][s['target']] in targets),None)
                chain.extend(['met:'+mid if mid else None,'rxn:'+step['reaction_id'] if step else None])
        side=max(2,math.ceil(len(chain)/4));entry=last_acetyl or [1,-3-margin]
        left,top=entry[0]+1,entry[1]+1
        ring=([(left,top-t) for t in range(side)]+[(left+t,top-side) for t in range(side)]+
              [(left+side,top-side+t) for t in range(side)]+[(left+side-t,top) for t in range(side)])
        for i,key in enumerate(chain):
            if key:
                put(key,ring[i*len(ring)//len(chain)],local)
        per_route['tca']=local;square=dict(left=left,top_offset=top,side_length=side,closed=route['closed'])
    fixed={k:v for k,v in offsets.items() if k in displayed_nodes}
    report_routes=[dict(r,node_offsets=per_route[r['id']],displayed_node_ids=[k for k in per_route[r['id']] if k in displayed_nodes],
                        displayed_reaction_ids=[rid for rid in r['reaction_ids'] if 'rxn:'+rid in displayed_nodes]) for r in ordered]
    return fixed,dict(routes=report_routes,square=square,color=RED,fixed_node_offsets=fixed,
                      glucose_top_margin=margin,geometry_scope='Displayed source-model nodes only; no missing edge is synthesized')
