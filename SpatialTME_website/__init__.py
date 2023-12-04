from flask import Flask
from flask_restful import Api
from celery_config import make_celery
from celery import Celery  # 确保导入 Celery


app = Flask(__name__)
app.config.from_object('eRNA_QTL.settings')
celery = Celery(app.import_name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)


api = Api(app)

app.url_map.strict_slashes = False

import eRNA_QTL.core
import eRNA_QTL.controllers
import eRNA_QTL.ajax
