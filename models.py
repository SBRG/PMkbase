from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class GrowthData(db.Model):
    """
    SQLAlchemy model for storing growth curve data for a given plate, species, and well.

    Attributes:
        id (int): Primary key.
        plateid (str): Plate identifier.
        specie (str): Species name.
        well (str): Well identifier.
        compound (str): Compound name.
        replicates (str): Number of replicates.
        signal_data (PickleType): List of signal values for the growth curve.
    """
    id = db.Column(db.Integer, primary_key=True)
    plateid = db.Column(db.String(50))
    specie = db.Column(db.String(50))
    well = db.Column(db.String(10))
    compound = db.Column(db.String(100))
    replicates = db.Column(db.String(10))
    signal_data = db.Column(db.PickleType)



class TraitData(db.Model):
    """
    SQLAlchemy model for storing trait (growth call) data for a given strain, plate, and well.

    Attributes:
        id (int): Primary key.
        strainid (str): Strain identifier.
        plateid (str): Plate identifier.
        specie (str): Species name.
        metadata_mods (str): Metadata or modifications.
        project (str): Project name.
        well (str): Well identifier.
        plate (str): Plate name.
        media (str): Media type.
        growth (int): Growth call (e.g., 0/1).
        compound (str): Compound name.
        desc (str): Description.
        strain (str): Strain name.
        phylo (str): Phylogroup or genome cluster.
        mlst (str): MLST type.
    """
    id = db.Column(db.Integer, primary_key=True)
    strainid = db.Column(db.String(50))
    plateid = db.Column(db.String(50))
    specie = db.Column(db.String(50))
    metadata_mods = db.Column(db.String(100))
    project = db.Column(db.String(50))
    well = db.Column(db.String(10))
    plate = db.Column(db.String(50))
    media = db.Column(db.String(50))
    growth = db.Column(db.Integer)
    compound = db.Column(db.String(100))
    desc = db.Column(db.String(100))
    strain = db.Column(db.String(50))
    phylo = db.Column(db.String(50))
    mlst = db.Column(db.String(50))



class KineticData(db.Model):
    """
    SQLAlchemy model for storing kinetic parameters for a given strain, plate, and well.

    Attributes:
        id (int): Primary key.
        plateid (str): Plate identifier.
        strainid (str): Strain identifier.
        strain (str): Strain name.
        specie (str): Species name.
        metadata_mods (str): Metadata or modifications.
        project (str): Project name.
        well (str): Well identifier.
        plate (str): Plate name.
        media (str): Media type.
        replicates (str): Number of replicates.
        compound (str): Compound name.
        keggid (str): KEGG compound ID.
        casid (str): CAS compound ID.
        maxresp (float): Maximum response value.
        maxresprate (float): Maximum response rate.
        timetill (float): Time until maximum response rate.
        auc (float): Area under the curve.
        growth (float): Growth value.
        mlst (str): MLST type.
        phylo (str): Phylogroup or genome cluster.
    """
    id = db.Column(db.Integer, primary_key=True)
    plateid = db.Column(db.String(50))
    strainid = db.Column(db.String(50))
    strain = db.Column(db.String(50))
    specie = db.Column(db.String(50))
    metadata_mods = db.Column(db.String(100))
    project = db.Column(db.String(50))
    well = db.Column(db.String(10))
    plate = db.Column(db.String(50))
    media = db.Column(db.String(50))
    replicates = db.Column(db.String(10))
    compound = db.Column(db.String(100))
    keggid = db.Column(db.String(100))
    casid = db.Column(db.String(100))
    maxresp = db.Column(db.Float)
    maxresprate = db.Column(db.Float)
    timetill = db.Column(db.Float)
    auc = db.Column(db.Float)
    growth = db.Column(db.Float)
    mlst = db.Column(db.String(50))
    phylo = db.Column(db.String(50))

