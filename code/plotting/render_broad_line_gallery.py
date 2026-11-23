"""Arrange existing profile and component-table images without regenerating them."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / 'output' / 'broad_lines'
EXISTING_TABLES = ROOT / 'output' / 'fit_diagnosis_20261014' / 'component_tables'
ADDITIONAL_TABLES = {
    'HC40NW': 'component_tables/HC40NW_components.png',
    'W44W': 'component_tables/W44W_components.png',
}
START_MARKER = '% BEGIN GENERATED BROAD-LINE GALLERY'
END_MARKER = '% END GENERATED BROAD-LINE GALLERY'


def image_panel(path):
    return '\n'.join([
        r'\begin{minipage}[t]{0.49\textwidth}',
        r'\vspace{0pt}\centering',
        r'\includegraphics[width=\linewidth]{' + path + '}',
        r'\end{minipage}',
    ])


def figure_row(paths, caption):
    panels = '\n\\hfill\n'.join(image_panel(path) for path in paths)
    return '\n'.join([
        r'\begin{figure}[!htbp]', r'\centering', panels,
        r'\caption{' + caption + '}', r'\end{figure}',
    ])


def build_gallery(profiles, paired_sources):
    lines = [START_MARKER, r'\clearpage', r'\normalsize',
             r'\section*{Molecular-line profiles and broad components}']
    for source in paired_sources:
        table_path = ADDITIONAL_TABLES.get(
            source, f'../fit_diagnosis_20261014/component_tables/{source}_components.png')
        paths = [f'figures/{source}.pdf', table_path]
        lines.append(figure_row(paths, source + ': molecular-line profiles (left) '
                                'and the broad-component table (right).'))
    remaining = sorted(set(profiles) - set(paired_sources))
    lines.extend([r'\clearpage', r'\subsection*{Additional molecular-line profiles}'])
    for index in range(0, len(remaining), 2):
        sources = remaining[index:index + 2]
        paths = [f'figures/{source}.pdf' for source in sources]
        lines.append(figure_row(paths, ' and '.join(sources) + ': molecular-line profiles.'))
    lines.extend([r'\clearpage', END_MARKER])
    return '\n'.join(lines) + '\n'


def main():
    profiles = sorted(path.stem for path in (OUTPUT / 'figures').glob('*.pdf'))
    paired_sources = [source for source in profiles
                      if source in ADDITIONAL_TABLES
                      or (EXISTING_TABLES / f'{source}_components.png').is_file()]
    gallery = build_gallery(profiles, paired_sources)
    document = OUTPUT / 'broad_line_table.tex'
    original = document.read_text()
    start = original.index(START_MARKER)
    end = original.index(END_MARKER, start) + len(END_MARKER)
    document.write_text(original[:start] + gallery + original[end:])
    print(f'Reused {len(paired_sources)} table/profile pairs and '
          f'{len(profiles) - len(paired_sources)} profile-only panels.')


if __name__ == '__main__':
    main()
