
from eRNA_QTL import app
from celery_config import make_celery

app.config.update(
    CELERY_BROKER_URL='redis://localhost:6379/0',
    CELERY_RESULT_BACKEND='redis://localhost:6379/0'
)

celery = make_celery(app)

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=1234)

    

