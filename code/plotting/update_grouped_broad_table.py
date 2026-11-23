"""Update the existing table, preserving the document's commented figure section."""
from pathlib import Path
import csv,json,math,re
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/broad_lines'
ORDER=['12CO10','12CO21','12CO43','HCO+10','HCO+32','HCN10','HCN32','SiO21']
LABEL={'12CO10':r'$^{12}$CO (1--0)','12CO21':r'$^{12}$CO (2--1)','12CO43':r'$^{12}$CO (4--3)','HCO+10':r'HCO$^+$ (1--0)','HCO+32':r'HCO$^+$ (3--2)','HCN10':'HCN (1--0)','HCN32':'HCN (3--2)','SiO21':'SiO (2--1)'}
def fmt(v,e):
 if e<=0:return f'${v:.3g}$'
 dec=min(5,max(0,1-math.floor(math.log10(e))))
 return f'${v:.{dec}f} \\pm {e:.{dec}f}$'
def render_table(raw):
 selected=[r for r in raw if float(r['amp'])>0 and float(r['fwhm'])>(7 if r['line']=='HCO+10' else 8)]
 sources=sorted({r['source'] for r in selected}); audit=[]
 heading=r'Target & Transition & $v_0$ & $T_p$ & FWHM & $\int T\,dv$ \\'
 units=r' & & (km s$^{-1}$) & (K) & (km s$^{-1}$) & (K km s$^{-1}$) \\'
 lines=[r'\setlength{\LTcapwidth}{\textwidth}',r'\begin{longtable}{cccccc}',r'\caption{Broad emission components grouped by HCO$^+$ (1--0) reference velocity, using the counterpart selection of Table~\ref{tab:broad-ratios}.}\label{tab:broad}\\',r'\toprule',heading,units,r'\midrule\endfirsthead',r'\multicolumn{6}{c}{Table \thetable\ continued}\\',r'\toprule',heading,units,r'\midrule\endhead',r'\midrule\multicolumn{6}{c}{Continued on next page}\\\endfoot',r'\bottomrule\endlastfoot']
 for si,source in enumerate(sources):
  entries=[r for r in selected if r['source']==source]
  anchors=sorted([r for r in entries if r['line']=='HCO+10'],key=lambda r:float(r['v']))
  groups=[]; used=set()
  for anchor in anchors:
   members=[anchor]
   for transition in [line for line in ORDER if line != 'HCO+10']:
    candidates=[r for r in entries if r['line']==transition and abs(float(r['v'])-float(anchor['v']))<=10]
    if candidates:
     members.append(min(candidates,key=lambda r:(abs(float(r['v'])-float(anchor['v'])),float(r['v']),int(r['component']))))
   groups.append(members)
   used.update((r['fit_file'],r['component']) for r in members)
  unmatched=[r for r in entries if (r['fit_file'],r['component']) not in used]
  blocks=[(f'Reference HCO$^+$ (1--0): $v_0={float(a["v"]):.2f}$ km s$^{{-1}}$',a,g) for a,g in zip(anchors,groups)]
  if unmatched:blocks.append(('Unmatched components' if anchors else 'No broad HCO$^+$ (1--0) reference',None,unmatched))
  if si:lines.append(r'\midrule[0.7pt]')
  for gi,(title,anchor,members) in enumerate(blocks):
   if gi:lines.append(r'\cmidrule[0.25pt](lr){2-6}')
   lines.append(r'\multicolumn{6}{c}{\footnotesize '+title+r'} \\*')
   members.sort(key=lambda r:(ORDER.index(r['line']),float(r['v']),int(r['component'])))
   for ri,r in enumerate(members):
    src=source+(r'$^{\dagger}$' if str(r['legacy']).lower()=='true' else '')
    cells=[src if ri==0 else '',LABEL[r['line']]]+[fmt(float(r[v]),float(r[e])) for v,e in [('v','ev'),('amp','eamp'),('fwhm','efwhm'),('integrated','eintegrated')]]
    last=gi==len(blocks)-1 and ri==len(members)-1
    lines.append(' & '.join(cells)+(r' \\' if last else r' \\*'))
    audit.append(dict(r,reference_component='' if anchor is None else anchor['component'],reference_v0='' if anchor is None else anchor['v'],velocity_difference='' if anchor is None else abs(float(r['v'])-float(anchor['v'])),group=gi+1,matched=anchor is not None))
 lines.extend([r'\end{longtable}',r'\noindent\footnotesize Notes. Only positive-emission components are included: HCO$^+$ (1--0) references require FWHM $>7$ km s$^{-1}$, and other components require FWHM $>8$ km s$^{-1}$. For each reference, the nearest qualifying component of each of $^{12}$CO (1--0), (2--1), (4--3), HCO$^+$ (3--2), HCN (1--0), (3--2), and SiO (2--1) is selected within $|v_0-v_{0,\mathrm{ref}}|\leq10$ km s$^{-1}$, exactly as in Table~\ref{tab:broad-ratios}; ties select the lower velocity, then component number. Matching is independent for each reference, so a counterpart may appear in more than one group. Components are not summed. All remaining broad components, including transitions outside this matching set, are listed separately. Targets without a qualifying HCO$^+$ reference remain ungrouped. Thick rules separate targets; thin rules separate groups. Errors are the reported fit uncertainties. $^{\dagger}$Legacy W44 products lack temperature-calibration metadata; values are preserved. Some other fits divide input labelled $T_{\rm mb}$ by a beam efficiency; absolute intensities require calibration verification. Velocity matching is a display convention, not an independent confirmation of physical association.'])
 return '\n'.join(lines)+'\n',audit

def main():
 raw=[]; versions={}
 for path in (ROOT/'output/line_fit').glob('*fit*.txt'):
  m=re.fullmatch(r'(.+)_(12CO\d+|HCN\d+|HCO[+p]\d+|SiO\d+)_fit(\d*)\.txt',path.name)
  if m:
   source,line,version=m.groups()
   source=source.upper() if source.lower() in {'w44w','w44n','w44c','w44s'} else source
   versions.setdefault((source,line.replace('HCOp','HCO+')),[]).append((int(version or 1),path.stat().st_mtime_ns,path))
 for (source,line),choices in versions.items():
  path=max(choices)[2]; text=path.read_text(); legacy='amp_K' not in text
  for row in text.splitlines():
   if not row.strip() or row.startswith('#'):continue
   a=list(map(float,row.split()))
   if legacy:
    v,ev,sigma,es,amp,ea,integ,ei=a[2:10]; width,ew=[math.sqrt(8*math.log(2))*x for x in (sigma,es)]
   else:amp,ea,v,ev,sigma,es,width,ew,integ,ei=a[2:12]
   if a[1]>0:
    raw.append(dict(source=source,line=line,fit_file=path.name,legacy=legacy,component=int(a[0]),amp=amp,eamp=ea,v=v,ev=ev,fwhm=width,efwhm=ew,integrated=integ,eintegrated=ei))
 fragment,audit=render_table(raw)
 document=OUT/'broad_line_table_standalone.tex';text=document.read_text()
 start=text.index(r'\setlength{\LTcapwidth}')
 end=text.index('% Figure section temporarily commented out;')
 preserved=text[end:]
 document.write_text(text[:start]+fragment+preserved)
 (OUT/'broad_line_table.tex').write_text(fragment)
 with (OUT/'grouped_emission_components.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(audit[0]));writer.writeheader();writer.writerows(audit)
 assert document.read_text().endswith(preserved)
 expected={(r['fit_file'],r['component']) for r in raw if r['amp']>0 and r['fwhm']>(7 if r['line']=='HCO+10' else 8)}
 assert {(r['fit_file'],r['component']) for r in audit}==expected
 ratios=json.loads((OUT/'broad_line_ratios.json').read_text())
 for record in ratios:
  actual={(r['line'],r['component']) for r in audit if r['source']==record['source'] and r['reference_component']==record['anchor']['component']}
  expected_group={(line,c['component']) for line,c in record['matches'].items() if c}
  compared={item for item in actual if item[0] in record['matches']}
  assert compared==expected_group,(record['source'],compared,expected_group)
 assert all(float(r['velocity_difference'])<=10 for r in audit if r['matched'])
 print(json.dumps({'rows':len(audit),'removed':len(raw)-len(audit),'matched_rows':sum(r['matched'] for r in audit),'unmatched_rows':sum(not r['matched'] for r in audit),'reference_groups':len({(r['source'],r['reference_component']) for r in audit if r['matched']})}))
if __name__=='__main__':main()
