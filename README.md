# AeroFEM

*Insert project description*

# Project Structure
```
├── aerofem
│   ├── main                        # Main logic for problem solving
│   │   └── example.py              # Example of a aerodynamic problem: planar static wing solved with Galerkin method. 
│   ├── models                      # Maing logic of models
│   │   └── aerodynamic            
│   │       └── lifting_line.py     
│   └── utils                       # Helper functions    
│       ├── aero_utils.py           
│       └── mesh_utils.py
├── README.md
└── setup.py
```

# Usage

- Clone repository
```bash
git clone git@github.com:glev1/aerofem.git
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
