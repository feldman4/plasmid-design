## Installation
```bash
pip install plasmid_design
```

## Example
In a new directory, save the following to `config.yaml`, then run `plasmid_design run`. The `drive_key` below comes from the URL of [this spreadsheet](https://docs.google.com/spreadsheets/d/1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM/edit#gid=52604569).

```yaml
drive_key: &drive_key 1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM
tables:
  # vector definitions
  templates: 
    table: drive_key:templates
    # selects just these rows
    gate: experiment == "example"
    drive_key: *drive_key

  # DNA and protein part definitions
  parts: 
    table: drive_key:parts
    drive_key: *drive_key
  
  # these DNA features will automatically be annotated
  features: 
    table: drive_key:features
    drive_key: *drive_key
  
  # these restriction sites will be avoided
  restriction_enzymes: 
    table: drive_key:enzymes
    # selects just these rows
    gate: pT02_BsaI == "x"
    drive_key: *drive_key

# optional key, override default arguments to DNAChisel
rt_args:
  species: h_sapiens # or e_coli, s_cerevisiae, etc
  seed: 1 # random seed
  k: 6 # avoid repeated DNA sequences of this length
```

The `run` command prints estimated cost and complexity scores. The following files will be generated:

```bash
gene_order.fa # genes, ready to order
gene_order_idt_scores.csv # IDT gene synthesis complexity scores (different from gblock complexity)
parts.csv # snapshot of input table
templates.csv # snapshot of input table
features.csv # snapshot of input table
reverse_translations/restriction_enzymes.csv # snapshot of input table
reverse_translations/input.fa # amino acid parts that need reverse translation
reverse_translations/output.fa # corresponding DNA, generated by DNA Chisel
reverse_translations/idt_scores.csv # per-part IDT gene synthesis complexity scores, not meaningful for short parts
reverse_translations/dna_chisel/ # DNA Chisel logs
```

## Configuration
Tables can be sourced in different ways.
- local csv
    - `table: /path/to/csv`
- private google sheets
    - `table: "drive:{spreadsheet name}/{sheet name}"`
- public google sheets
    - `table: "drive_key:{sheet_name}"`
    - `drive_key: {drive_key}`

For the public option, first enable link access in google sheets, then copy `drive_key` from the URL: `https://docs.google.com/spreadsheets/d/{drive_key}/`.

## Changing the vector backbone
The DNA sequences upstream and downstream of the designed gene in the plasmid map are identified by having "up" or "down" in the part name. To change the vector backbone, define new upstream and downstream DNA parts, then change the template column in the templates sheet accordingly.

