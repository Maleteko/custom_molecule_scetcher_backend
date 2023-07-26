# Django Backend for Molecular 2D Structure Visualization

This Django backend is part of a project that provides a web-based tool for visualizing the 2D structure of molecules based on SMILES strings. The backend is responsible for handling incoming SMILES data, calculating the 2D coordinates of the molecules using RDKit, and responding with the calculated coordinates in a JSON format.

## Prerequisites

Before running the backend, make sure you have the following installed:

- Python 3.x
- Django
- RDKit

## Installation

1. Clone this repository to your local machine.
2. Navigate to the project's root directory.

```bash
cd path/to/django-backend/
```

3. Create a virtual environment (optional but recommended).

```bash
python -m venv venv
source venv/bin/activate
```

4. Install the required Python packages.

```bash
pip install -r requirements.txt
```

## Configuration

1. Open the `settings.py` file located in the `molecule_visualizer` directory.
2. Configure the database settings if required (default is SQLite).
3. Ensure that the `ALLOWED_HOSTS` setting includes the host from which the Vue.js frontend will send API requests.

## Running the Server

To start the Django development server, run the following command:

```bash
python manage.py runserver
```

The backend server will now be running at `http://localhost:8000/`.

## API Endpoint

The backend exposes a single API endpoint for handling SMILES data:

- **URL:** `/coordinates/`
- **Method:** POST
- **Request Body:**

```json
{
  "smiles": "YOUR_SMILES_STRING"
}
```

- **Response:**

```json
{
  "atoms": [
    {"element": "ELEMENT_1", "x": X_COORD_1, "y": Y_COORD_1},
    {"element": "ELEMENT_2", "x": X_COORD_2, "y": Y_COORD_2},
    ...
  ],
  "bonds": [
    {"atom1": ATOM_INDEX_1, "atom2": ATOM_INDEX_2, "order": BOND_ORDER},
    ...
  ]
}
```

## Dependencies

The backend utilizes the following main dependencies:

- Django: Web framework for building the backend.
- RDKit: Chemistry toolkit for calculating 2D coordinates of molecules from SMILES strings.

## Contributing

Contributions to this project are welcome. If you find any issues or want to suggest improvements, feel free to open an issue or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).
