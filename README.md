## Installation
`pip install plasmid_design`

## Example
In a new directory, save the following to `config.yaml`, then run `plasmid_design run`. The configuration below refers to [this spreadsheet](https://docs.google.com/spreadsheets/d/1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM/edit#gid=52604569).

```yaml
tables:
  # vector definitions
  templates: 
    table: drive_key:templates
    # selects just these rows
    gate: experiment == "example"
    drive_key: 1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM

  # DNA and protein part definitions
  parts: 
    table: drive_key:parts
    drive_key: 1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM
  
  # these DNA features will automatically be annotated
  features: 
    table: drive_key:features
    drive_key: 1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM
  
  # these restriction sites will be avoided
  restriction_enzymes: 
    table: drive_key:enzymes
    # selects just these rows
    gate: pT02_BsaI == "x"
    drive_key: 1QWFQUlIJYERJ6zY-THD9uNagD2a_He7RUfcOeF2eEKM
```

The `run` command prints estimated cost and complexity scores. the following files will be generated.

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
The DNA sequences upstream and downstream of the designed gene are identified by having "up" or "down" in the part name. To change the vector backbone, define new upstream and downstream DNA parts, then change the template column in the templates sheet accordingly.

