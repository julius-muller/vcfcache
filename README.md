

``` bash
vepstash/
├── config/
│   └── vep_defaults.yaml
├── docker/
│   ├── Dockerfile              # Base image with VEP, Python, Nextflow
│   ├── docker-compose.yml      # Service definition
│   └── requirements.txt        # Python dependencies
├── vepstash/
│   ├── __init__.py
│   ├── cache.py       # database creation/management
│   ├── convert.py     # annotation using cache
│   └── config.py      # yaml config handling
├── workflow/
│   ├── main.nf        # nextflow pipeline
│   ├── nextflow.config
│   └── modules/       # nextflow process modules
│       ├── annotate.nf
│       ├── intersect.nf
│       └── vep.nf
├── tests/
│   ├── test_cache.py
│   ├── test_convert.py
│   └── data/         # test fixtures
├── scripts/
│   └── entrypoint.sh  # docker entrypoint
├── README.md
├── pyproject.toml
├── setup.py
└── .dockerignore
```