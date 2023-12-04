# celery_config.py
from celery import Celery

def make_celery(app):
    # 实例化一个新的 Celery 对象
    celery = Celery(app.import_name, broker=app.config['CELERY_BROKER_URL'])

    # 更新 Celery 配置
    celery.conf.update(app.config)

    # (可选) 在 Flask 应用上下文中运行 Celery 任务
    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    return celery

