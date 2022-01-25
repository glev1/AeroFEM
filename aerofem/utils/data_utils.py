import os
import json
import pickle


def set_main_path() -> str:
    """Find main path to the project"""

    return os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.abspath(__file__))))

def get_param(config_name: str, value_name: str):
    """Get parameters values.
    Args:
        config_name (str): Name of config file.
        config_key (str): Name of the key within the
        parameter being loaded.
    Returns:
        The value corresponding to the requested key.
    """
    main_path = set_main_path()

    file = os.path.join(main_path, f"aerofem/conf/{config_name}.json")

    with open(file, encoding='utf-8') as f:
        data = json.load(f)
    return data[value_name]

def save_object(obj, filename):

    main_path = set_main_path()
    file = os.path.join(main_path, f"projects/{filename}.pkl")

    with open(file, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def load_object(path, filename):

    main_path = set_main_path()
    file = os.path.join(main_path, f"{path}/{filename}.pkl")

    with open(file, 'rb') as inp:  # Overwrites any existing file.
        obj = pickle.load(inp)

    return obj