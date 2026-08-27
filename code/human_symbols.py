#!/usr/bin/env python3
"""
Shared human-symbol -> GSE326743-panel-symbol resolver.

Extracted from build_genesets_human.py (Phase 7 Job A) so that the Job-A tier build and the
spleen marker build resolve symbols identically. Behaviour is unchanged; see the note below on
why the naive alias union is wrong.

Two authorities, both pinned files, never memory:
  genesets/h1_candidate/GSE326743_gene_panel_5093.csv   the panel (gene_name, gene_id)
  genesets/hgnc_pin/hgnc_symbol_alias_ensembl.*.csv.gz  HGNC symbol/alias/prev -> Ensembl

Folding approved symbols, previous symbols and aliases into ONE dict silently produces wrong
genes -- CCNL1 collapses onto MCM2, MIF onto AMH, OPA1 onto MED12 -- because those names appear
in some other gene's alias list. The rule that prevents it: A SYMBOL THAT IS ITSELF AN APPROVED
HGNC SYMBOL ONLY EVER RESOLVES TO ITS OWN ENSEMBL ID. Aliases are consulted only for names HGNC
does not recognise as approved, and only when the alias is unambiguous.
"""
import csv, gzip

GS = '/workspace/genesets'
PANEL_CSV = GS + '/h1_candidate/GSE326743_gene_panel_5093.csv'
HGNC_CSV  = GS + '/hgnc_pin/hgnc_symbol_alias_ensembl.2026-08-27.csv.gz'

_PANEL_ROWS = list(csv.DictReader(open(PANEL_CSV)))
PANEL    = {r['gene_name'] for r in _PANEL_ROWS}
ENSG2SYM = {r['gene_id']: r['gene_name'] for r in _PANEL_ROWS}

_HG = list(csv.DictReader(gzip.open(HGNC_CSV, 'rt')))
APPROVED, PREV, ALIAS = {}, {}, {}
_prev_multi, _alias_multi = set(), set()
for _r in _HG:
    if _r['status'] != 'Approved' or not _r['ensembl_gene_id']:
        continue
    APPROVED[_r['symbol']] = _r['ensembl_gene_id']
for _r in _HG:
    if _r['status'] != 'Approved' or not _r['ensembl_gene_id']:
        continue
    _e = _r['ensembl_gene_id']
    for _n in [x.strip() for x in (_r['prev_symbol'] or '').replace('"', '').split('|') if x.strip()]:
        if PREV.setdefault(_n, _e) != _e:
            _prev_multi.add(_n)
    for _n in [x.strip() for x in (_r['alias_symbol'] or '').replace('"', '').split('|') if x.strip()]:
        if ALIAS.setdefault(_n, _e) != _e:
            _alias_multi.add(_n)
for _n in _prev_multi:
    PREV.pop(_n, None)
for _n in _alias_multi:
    ALIAS.pop(_n, None)

RESOLUTIONS = []          # (asked, resolved_to, route) -- appended as a side effect of resolve()


def resolve(sym):
    """Human symbol -> the symbol the GSE326743 panel actually uses, or None if off-panel."""
    if sym in PANEL:
        return sym
    if sym in APPROVED:                        # approved symbol: ONLY its own Ensembl id counts
        e = APPROVED[sym]
        if e in ENSG2SYM:
            RESOLUTIONS.append((sym, ENSG2SYM[e], 'HGNC approved symbol -> %s -> panel legacy symbol' % e))
            return ENSG2SYM[e]
        return None
    for tag, idx in (('HGNC previous symbol', PREV), ('HGNC unambiguous alias', ALIAS)):
        e = idx.get(sym)
        if e and e in ENSG2SYM:
            RESOLUTIONS.append((sym, ENSG2SYM[e], '%s -> %s' % (tag, e)))
            return ENSG2SYM[e]
    return None


def known(sym):
    """True if HGNC recognises the string at all (approved, previous, or unambiguous alias)."""
    return sym in PANEL or sym in APPROVED or sym in PREV or sym in ALIAS


def on_panel(genes):
    return sorted({r for r in (resolve(g) for g in genes) if r})
