# AeroFEM

*Insert project description*

# Project Structure
```
├── aerofem
│   ├── conf                       # Model and parameter configuration
│   │   ├── example_mesh.json
│   │   └── example_wing.json
│   ├── __init__.py
│   ├── models                     # Main code of models
│   │   ├── aerodynamic
│   │   │   ├── __init__.py
│   │   │   └── lifting_line.py
│   │   └── __init__.py
│   └── utils                      # Helper functions
│       ├── aero_utils.py
│       ├── data_utils.py
│       ├── __init__.py
│       └── mesh_utils.py
├── Makefile                   
├── projects                       # Saved projects
│   └── example.pkl
├── README.md
├── requirements.txt               # Dependencies
└── tutorials                      # Usage examples
    └── example.py
```

# Usage

- Clone repository
```bash
git@github.com:glev1/AeroFEM.git
```

- Create virtualenv
```bash
python3 -m venv ~/.aerofem && source ~/.aerofem/bin/activate
```

- Install dependencies 
```bash
make install 
```

Run simple example:
```bash
python3 -m tutorials.example
```


# References

*Insert articles*
