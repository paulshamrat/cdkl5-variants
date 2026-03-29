# 03 Binding

`03_binding` is the planned modular binding stage of TGVR.

## Intended Scope

This stage is meant to:

- prepare partner-ready mutation input files
- ingest returned external binding predictor outputs
- summarize partner-specific binding DDG values

It is conceptually downstream of `01_variant_curation`.

## Important Current Reality

The older stage notes describe a `run_binding.py` entrypoint and a binding implementation layer, but those public files are not currently present under [workflow/tgvr/scripts](/home/paul/cdkl5-variants/workflow/tgvr/scripts) or [workflow/tgvr/src](/home/paul/cdkl5-variants/workflow/tgvr/src).

So this stage should currently be treated as:

- documented design intent
- archived workflow direction
- not yet a public runnable TGVR stage in the present tree

## Expected Future Responsibilities

When completed, this stage is expected to cover:

- input preparation for predictors such as `SAAMBE-3D`, `FoldX`, `mCSM-PPI2`, and `DDMutPPI`
- ingestion of returned external result files
- binding-DDG summaries for selected label groups

## Relationship To Legacy Work

This stage continues the legacy binding idea in a more modular form, but the current repository state has not yet exposed the runnable public binding interface alongside the other TGVR scripts.
