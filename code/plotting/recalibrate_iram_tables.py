"""Update saved table selections and ratios without changing component matching."""
from pathlib import Path
import csv
import json
import math
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'code/other'))
from astropy.io import fits
from recalibrate_iram_products import VERSION, metadata, resolveInput
from beam_efficiency import iram_header_efficiency
from rebuild_iram_products import OUT, backup, scaleWindow, drawWindow
import build_broad_line_products as profiles

PAIRS = [('12CO21','12CO43'), ('HCO+10','12CO43'), ('HCO+10','12CO21'), ('HCO+10','12CO10')]
NAMES = ['C21/C43','HCO+10/C43','HCO+10/C21','HCO+10/C10']
SCALE_KEYS = {'amp','eamp','integrated','eintegrated','rms','integrated_K_kms',
              'error_K_kms','rms_K','measurement_error_K_kms','calibration_error_K_kms'}


def factorFor(source, line, fitFile=None):
    if 'W44' in source.upper():
        return 1.0
    if fitFile:
        path = ROOT / fitFile
        if not path.exists():
            path = ROOT/'code/fitting/line_fit'/Path(fitFile).name
        if path.exists():
            meta = metadata(path.read_text())
            if 'input_file' in meta:
                inputPath = resolveInput(meta['input_file'])
                return iram_header_efficiency(fits.getheader(inputPath), inputPath)
    path = ROOT/'snr_data'/f'{source}_{line}.fits'
    return iram_header_efficiency(fits.getheader(path), path) if path.exists() else 1.0


def scaleNested(value, factor):
    if isinstance(value, dict):
        for key, item in value.items():
            if key in SCALE_KEYS and isinstance(item,(float,int)):
                value[key] *= factor
            elif isinstance(item,(dict,list)):
                scaleNested(item, factor)
    elif isinstance(value,list):
        for item in value:
            scaleNested(item,factor)


def ratio(a,b):
    if not a or not b:
        return None
    x,dx,y,dy = a.get('integrated_K_kms'),a.get('error_K_kms'),b.get('integrated_K_kms'),b.get('error_K_kms')
    if any(v is None or not math.isfinite(v) for v in (x,dx,y,dy)) or y == 0:
        return None
    return dict(value=x/y,error=math.hypot(dx/y,x*dy/y**2))


def updateRatioTables():
    path = ROOT/'output/line_ratio_tables/line_ratios_audit.json'
    doc = json.loads(path.read_text())
    for method,records in doc.items():
        for record in records:
            if record.get('calibration_revision') != VERSION:
                if method == 'window':
                    scaleWindow(record)
                else:
                    for line,entry in record['lines'].items():
                        scaleNested(entry,factorFor(record['source'],line,entry.get('fit_file')))
                record['calibration_revision'] = VERSION
            record['ratios'] = {name:ratio(record['lines'].get(a),record['lines'].get(b))
                                for name,(a,b) in zip(NAMES,PAIRS)}
        texPath = path.parent/f'{method}_ratios.tex'
        text = texPath.read_text()
        lines=[]
        for record in records:
            values = ['-' if record['ratios'][name] is None else
                       profiles.fmt(record['ratios'][name]['value'],record['ratios'][name]['error']) for name in NAMES]
            lines.append(' & '.join([record['source']]+values+['','',''])+r' \\')
        text = text.split(r'\startdata')[0]+r'\startdata'+'\n'+'\n'.join(lines)+'\n'+r'\enddata'+text.split(r'\enddata')[1]
        backup(texPath);texPath.write_text(text)
    backup(path);path.write_text(json.dumps(doc,indent=2)+'\n')
    csvPath=path.parent/'line_ratios.csv';backup(csvPath)
    fields=['method','source']+[key for n in NAMES for key in (n,n+'_error')]+['FWHM','SNR age','H2 luminosity']
    with csvPath.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for method,records in doc.items():
            for record in records:
                row=dict(method=method,source=record['source'])
                for name,item in record['ratios'].items():
                    row[name]='-' if item is None else item['value']
                    row[name+'_error']='-' if item is None else item['error']
                writer.writerow(row)
    # Keep the table's historical window selections in their own directory.
    directory=ROOT/'output/line_ratio_tables/window_inputs'
    for record in doc['window']:
        target=directory/f"{record['source']}_window.json"
        backup(target);target.write_text(json.dumps(record,indent=2)+'\n')
        drawWindow(record,directory)


def updateComponentCsv():
    paths = list((ROOT/'output/broad_lines/single_source').glob('*_components.csv'))
    paths += [ROOT/'output/broad_lines/grouped_emission_components.csv']
    for path in paths:
        if 'W44' in path.name.upper():
            continue
        backupPath=OUT/'before'/path.relative_to(ROOT)
        sourcePath=backupPath if backupPath.exists() else path
        rows=list(csv.DictReader(sourcePath.open()))
        for row in rows:
            factor=factorFor(row['source'],row['line'],row.get('fit_file'))
            if factor != 1:
                for key in SCALE_KEYS & row.keys():
                    if row[key]: row[key]=f'{float(row[key])*factor:.12g}'
        backup(path)
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def updateTex():
    paths=[ROOT/'output/broad_lines/broad_line_table.tex',ROOT/'output/broad_line_table.tex']
    paths+=list((ROOT/'output/broad_lines/component_tables').glob('*.tex'))
    paths+=list((ROOT/'output/fit_diagnosis_20261014/component_tables').glob('*.tex'))
    for path in paths:
        if 'W44' in path.name.upper():continue
        text=path.read_text()
        if VERSION in text:continue
        source=path.name.split('_components')[0] if '_components' in path.name else ''
        lines=[]
        for line in text.splitlines():
            cells=line.split(' & ')
            if len(cells)==6:
                first=re.sub(r'\\qflag\{[^}]*\}|\$.*?\$', '',cells[0]).strip()
                if first and not first.startswith('\\') and first not in ('Target','Source'):
                    source=first
                transition=re.sub(r'\\qflag\{[^}]*\}','',cells[1])
                # Scientific labels retain molecule and upper/lower transition digits.
                compact=re.sub(r'[^A-Za-z0-9+]', '',transition).lower()
                mol='HCO+' if 'hco' in compact else 'HCN' if 'hcn' in compact else 'SiO' if 'sio' in compact else '12CO' if 'co' in compact else ''
                digits=re.findall(r'\d',compact)
                key=mol+''.join(digits[-2:]) if len(digits)>=2 else ''
                factor=factorFor(source,key) if key else 1
                if factor!=1:
                    for index in (3,5):
                        match=re.search(r'\$\s*([-+\d.eE]+)\s*\\pm\s*([-+\d.eE]+)\s*\$',cells[index])
                        if match:
                            cells[index]=cells[index].replace(match[0],profiles.fmt(float(match[1])*factor,float(match[2])*factor))
                    line=' & '.join(cells)
            lines.append(line)
        backup(path);path.write_text('% '+VERSION+'; IRAM only, W44 retained\n'+'\n'.join(lines)+'\n')


def updateBroadRatios():
    path=ROOT/'output/broad_lines/broad_line_ratios.json'
    records=json.loads(path.read_text())
    pairs=[('HCO+10','12CO10'),('HCO+10','12CO21'),('12CO21','12CO10'),
           ('12CO43','12CO10'),('12CO43','12CO21'),('HCO+10','HCN10'),('HCO+10','12CO43')]
    for record in records:
        if record.get('calibration_revision')==VERSION:continue
        for line,c in record['matches'].items():
            if c:
                factor=factorFor(record['source'],line,c.get('fit_file'))
                for key in ('integrated','error'):c[key]*=factor
        record['anchor']=record['matches']['HCO+10'].copy()
        record['ratios']=[record['matches'][a]['integrated']/record['matches'][b]['integrated']
                          if record['matches'][a] and record['matches'][b] and record['matches'][b]['integrated']>0 else None for a,b in pairs]
        record['calibration_revision']=VERSION
    backup(path);path.write_text(json.dumps(records,indent=2)+'\n')
    for target in (path.with_suffix('.tex'),ROOT/'output/broad_lines/broad_line_table.tex'):
        text=target.read_text();lines=[]
        for line in text.splitlines():
            cells=line.split(' & ')
            if len(cells)==10:
                source=cells[0].strip()
                match=next((r for r in records if r['source']==source and cells[1].strip()==f"{r['anchor']['v']:.2f}"),None)
                if match:
                    cells[3:]=['-' if v is None else f'{v:.3f}' for v in match['ratios']]
                    line=' & '.join(cells)+r' \\'
            lines.append(line)
        backup(target);target.write_text('\n'.join(lines)+'\n')


def main():
    updateRatioTables()
    updateComponentCsv()
    updateTex()
    updateBroadRatios()
    print('Updated IRAM table values and ratios; historical selections preserved.')


if __name__=='__main__':
    main()
