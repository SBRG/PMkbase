from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class GrowthData(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    plateid = db.Column(db.String(50))
    specie = db.Column(db.String(50))
    well = db.Column(db.String(10))
    compound = db.Column(db.String(100))
    replicates = db.Column(db.String(10))
    signal_data = db.Column(db.PickleType)



class TraitData(db.Model):
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

