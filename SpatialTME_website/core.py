from flask_pymongo import PyMongo

from eRNA_QTL import app

app.config['MONGO_URI']="mongodb://0.0.0.0:27017/eRNA_QTL"
mongo = PyMongo(app)