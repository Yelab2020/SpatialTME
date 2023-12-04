from __future__ import division
import flask_restful
from flask import Flask, send_from_directory
import pymongo
from pymongo.collation import Collation
from eRNA_QTL import app, api, celery
from eRNA_QTL.core import mongo
from flask_restful import Resource, fields, marshal_with, reqparse, marshal,request
from flask import url_for,redirect,current_app
from collections import Counter, defaultdict
import werkzeug
########
import subprocess,string,tempfile,shlex,multiprocessing,uuid,json,shutil,os,re,zipfile
from pprint import pprint
####
from scipy import stats
import numpy as np
import math
import subprocess
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import tarfile
import zipfile
import shutil
import pandas as pd
import rpy2.robjects as robjects
import re
import logging


def command_excute(command):
    args = shlex.split(command)
    #subprocess.check_call(args)
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(args, stdout=devnull, stderr=devnull)

### this function for key not in dict
def get_values(x,y):
    if x in y.keys():
        return y[x]
    else:
        return ''

###########
# Samples   ###
###########
samples_term ={
    "sample_id": fields.String,
    "dataset_name": fields.String,
    "cancer_type": fields.String,
    "organ": fields.String,
    "organ_detail": fields.String,
    "sample_status": fields.String,
    "boundary_define": fields.String,
    "published_date": fields.String,
    "accessible_ID": fields.String,
    "platform": fields.String,
    "doi": fields.String,
    "orginal_description": fields.String,
}
Sample_list = {
    'Sample_list':fields.List(fields.Nested(samples_term))
}

class Sample_search(Resource):
    @marshal_with(Sample_list)
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('sample_id_flask') #接受参数
        parser.add_argument('dataset_name_flask')
        parser.add_argument('cancer_type_flask')
        parser.add_argument('organ_flask')
        parser.add_argument('organ_detail_flask')
        parser.add_argument('sample_status_flask')
        parser.add_argument('boundary_define_flask')
        parser.add_argument('published_date_flask')
        parser.add_argument('accessible_ID_flask')
        parser.add_argument('platform')
        parser.add_argument('doi_flask')
        parser.add_argument('orginal_description_flask')
        args = parser.parse_args()
        condition = []
        if args['cancer_type_flask']!='':
            condition.append({'$match': {'cancer_type': args['cancer_type_flask']}})
        if args['organ_flask']!='':
            condition.append({'$match': {'organ': args['organ_flask']}})
        if args['accessible_ID_flask']!='':
            condition.append({'$match': {'accessible_ID': args['accessible_ID_flask']}})
        if args['platform']!='':
            condition.append({'$match': {'platform': args['platform']}})    
        if args['sample_id_flask']!='':
            condition.append({'$match': {'sample_id': args['sample_id_flask']}})
        print(condition)
        res = list(mongo.db.samples.aggregate(condition))
        return {'Sample_list': res}
api.add_resource(Sample_search, '/api/samples')

###########
# Download demo  ###
###########
class DownloadFile_scRNA_mat(Resource):
    def get(self):
        try:
            directory = "/work/sjt/SpatialTME/SpatialTME_localVersion/eRNA_QTL/static/eRNA_QTL/demo_data"  # 文件所在目录
            filename = "scRNAref_mat_demo.rds.gz"
            return send_from_directory(directory, filename, as_attachment=True)
        except Exception as e:
            # 这里可以更详细地处理不同类型的异常
            return {"error": str(e)}, 500
api.add_resource(DownloadFile_scRNA_mat, '/download_scrna_mat')

class DownloadFile_scRNA_metadata(Resource):
    def get(self):
        try:
            directory = "/work/sjt/SpatialTME/SpatialTME_localVersion/eRNA_QTL/static/eRNA_QTL/demo_data"  # 文件所在目录
            filename = "scRNAref_metadata_demo.csv"
            return send_from_directory(directory, filename, as_attachment=True)
        except Exception as e:
            # 这里可以更详细地处理不同类型的异常
            return {"error": str(e)}, 500
api.add_resource(DownloadFile_scRNA_metadata, '/download_scrna_metadata')

class DownloadFile_ST(Resource):
    def get(self):
        try:
            directory = "/work/sjt/SpatialTME/SpatialTME_localVersion/eRNA_QTL/static/eRNA_QTL/demo_data"  # 文件所在目录
            filename = "TumorST_demo.tar.gz"
            return send_from_directory(directory, filename, as_attachment=True)
        except Exception as e:
            # 这里可以更详细地处理不同类型的异常
            return {"error": str(e)}, 500
api.add_resource(DownloadFile_ST, '/download_st')

###########
# Upload   ###
###########
BASE_UPLOAD_FOLDER = '/work/sjt/SpatialTME/SpatialTME_localVersion/eRNA_QTL/static/upload'



class FileUpload(Resource):
    def post(self):
        try:
            if 'stFile' not in request.files or 'scrnaMetadataFile' not in request.files or 'scrnaMatFile' not in request.files:
                return {'message': 'Please fill out all fields.'}, 400
            stFile = request.files['stFile']
            scrnaMatFile = request.files['scrnaMatFile']
            scrnaMetadataFile = request.files['scrnaMetadataFile']
            tumorCell = request.form['tumorCell']
            epithelialCell = request.form['epithelialCell']
            stromalCell = request.form['stromalCell']
            stromal_cells = re.split(r'[,\s/]+', stromalCell)
            print(f"stromal_cells: {stromal_cells}")
            email = request.form['email']
            print(f"email: {email}")

            email_folder = os.path.join(BASE_UPLOAD_FOLDER, werkzeug.utils.secure_filename(email))
            print(f"email_folder: {email_folder}")

            if os.path.exists(email_folder):
                return {'message': 'A process with the same email is already running. Please wait for it to complete before starting a new one.'}, 400
            os.makedirs(email_folder, exist_ok=True)

            st_filename = werkzeug.utils.secure_filename(stFile.filename)
            st_file_path = os.path.join(email_folder, st_filename)
            stFile.save(st_file_path)

            scrna_mat_filename = werkzeug.utils.secure_filename(scrnaMatFile.filename)
            scrna_mat_file_path = os.path.join(email_folder, scrna_mat_filename)
            scrnaMatFile.save(scrna_mat_file_path)

            scrna_metadata_filename = werkzeug.utils.secure_filename(scrnaMetadataFile.filename)
            scrna_metadata_file_path = os.path.join(email_folder, scrna_metadata_filename)
            scrnaMetadataFile.save(scrna_metadata_file_path)

            # 解压 ST 文件
            extract_success, extract_message = self.extract_files(st_file_path, email_folder)
            if not extract_success:
                raise Exception(f'ST file extraction failed: {extract_message}')

            # 验证ST文件是否符合要求
            validate_success, validate_message = self.validate_st_contents(email_folder)
            if not validate_success:
                raise Exception(f'ST file content does not meet requirements: {validate_message}')

            # 检查 scrna 文件是否符合要求
            validate_success_scrna, validate_message_scrna = self.validate_scrna_files(scrna_mat_file_path, scrna_metadata_file_path, tumorCell, epithelialCell, stromal_cells)
            if not validate_success_scrna:
                raise Exception(validate_message_scrna)

            return {'message': 'File upload successful'}

        except Exception as e:
            # 如果发生异常，删除 email_folder
            if os.path.exists(email_folder):
                shutil.rmtree(email_folder)
            return {'message': f'File upload failed: {str(e)}'}, 400


    def extract_files(self, file_path, dest_folder):
        try:
            if file_path.endswith('.tar.gz') or file_path.endswith('.tar'):
                with tarfile.open(file_path, 'r:*') as tar:
                    tar.extractall(path=dest_folder)
            elif file_path.endswith('.zip'):
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    zip_ref.extractall(dest_folder)
            else:
                logging.error(f"Unsupported file format: {file_path}")
                return False, "Unsupported file format"
            return True, "Success"
        except Exception as e:
            logging.error(f"Error extracting file {file_path}: {e}")
            return False, str(e)


    def validate_st_contents(self, folder):
        # 初始化标志
        matrix_exists = False
        spatial_exists = False

        # 遍历目录树来搜索需要的文件和文件夹
        for subdir, dirs, files in os.walk(folder):
            for file in files:
                # 检查 matrix 文件或文件夹
                if file == 'filtered_feature_bc_matrix.h5':
                    matrix_path = os.path.join(subdir, file)
                    matrix_exists = True
                # 检查 spatial 文件夹中必需的文件
                if 'spatial' in dirs:
                    spatial_folder = os.path.join(subdir, 'spatial')
                    image_files_exist = any(os.path.exists(os.path.join(spatial_folder, f)) for f in ['tissue_hires_image.png', 'tissue_lowres_image.png'])
                    scale_factors_exist = os.path.exists(os.path.join(spatial_folder, 'scalefactors_json.json'))
                    positions_list_exist = os.path.exists(os.path.join(spatial_folder, 'tissue_positions_list.csv'))
                    if image_files_exist and scale_factors_exist and positions_list_exist:
                        spatial_exists = True
                if matrix_exists and spatial_exists:
                    break  # 如果已找到所需的所有文件，则无需进一步搜索

            if matrix_exists and spatial_exists:
                break  # 如果已找到所需的所有文件，则无需进一步搜索

        # 组合缺失的文件信息
        missing_files = []
        if not matrix_exists:
            missing_files.append("filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix directory")
        if not spatial_exists:
            missing_files.append("required files in spatial directory")

        # 如果缺失文件，返回 False 和缺失信息
        if missing_files:
            return False, f"Missing files: {', '.join(missing_files)}"

        # 如果一切正常，返回 True
        return True, "Success"


    def validate_scrna_files(self, scrna_mat_file_path, scrna_metadata_file_path, tumorCell, epithelialCell, stromal_cells):
    # 初始化DataFrame为空或者其他默认值
        #scrna_mat_df = pd.DataFrame()  # 或者适当的初始值
        #scrna_metadata_df = pd.DataFrame()  # 或者适当的初始值

        try:
            # 读取 scrnaMetadataFile
            scrna_metadata_df = pd.read_csv(scrna_metadata_file_path, index_col=0)
            if scrna_mat_file_path.endswith('.rds') or scrna_mat_file_path.endswith('.rds.gz'):
                readRDS = robjects.r['readRDS']
                scrna_mat_df = readRDS(scrna_mat_file_path)
                colnames = robjects.r['colnames']
                pass
            elif scrna_mat_file_path.endswith('.csv') or scrna_mat_file_path.endswith('.csv.gz'):
                scrna_mat_df = pd.read_csv(scrna_mat_file_path, index_col=0)
            # 如果没有合适的分支来读取文件，应该抛出一个错误
            else:
                raise ValueError("Unsupported file format for scRNA matrix file")

            # 检查行名和列名是否一致
            if set(scrna_metadata_df.index).issubset(set(scrna_mat_df.colnames)) and set(scrna_mat_df.colnames).issubset(set(scrna_metadata_df.index)):
                # 检查第一列的值是否包含 tumorCell、epithelialCell 和 stromalCell
                cell_types = [tumorCell, epithelialCell] + stromal_cells
                if set(cell_types).issubset(set(scrna_metadata_df.iloc[:, 0])):
                    return True, "scRNA matrix file and metadata file are valid."
                else:
                    return False, "scRNA metadata file does not contain specified cell types."
            else:
                return False, "scRNA matrix file column names do not match the row names in scRNA metadata file."

        except Exception as e:
            # 如果读取文件失败，处理异常
            return False, str(e)
    


 

api.add_resource(FileUpload, '/FileUpload')

###########
# Upload-check   ###
###########



###########
# Run   ###
###########
class RunAnalysis(Resource):
    def post(self):
        stFile = request.files['stFile']
        scrnaMatFile = request.files['scrnaMatFile']
        scrnaMetadataFile = request.files['scrnaMetadataFile']
        tumorCell = request.form['tumorCell']
        epithelialCell = request.form['epithelialCell']
        stromalCell = request.form['stromalCell']
        email = request.form['email']

        email_folder = os.path.join(BASE_UPLOAD_FOLDER, werkzeug.utils.secure_filename(email))
        if not os.path.exists(email_folder):
            return {'message': '未找到相应的上传文件夹'}, 404
        
        st_filename = werkzeug.utils.secure_filename(stFile.filename)
        scrna_mat_filename = werkzeug.utils.secure_filename(scrnaMatFile.filename)
        scrna_metadata_filename = werkzeug.utils.secure_filename(scrnaMetadataFile.filename)

        #  run_r_script 函数调用 R 脚本，并传入文件夹路径
        #result = run_r_script(email_folder, st_filename, scrna_mat_filename, scrna_metadata_filename, tumorCell,epithelialCell,stromalCell,email)
        task = run_r_script_task.delay(email_folder, st_filename, scrna_mat_filename, scrna_metadata_filename, tumorCell, epithelialCell, stromalCell, email)
        return {'task_id': task.id}, 202

@celery.task        
def run_r_script_task(email_folder, st_filename, scrna_mat_filename, scrna_metadata_filename, tumorCell, epithelialCell, stromalCell,email):
    command = ['Rscript', '/work/sjt/SpatialTME/ST_pipeline/PipeLine.R', str(email_folder), str(st_filename),str(scrna_mat_filename),str(scrna_metadata_filename),str(tumorCell),str(epithelialCell),str(stromalCell),str(email)]
    result = subprocess.run(command)
    return result.stdout

api.add_resource(RunAnalysis, '/RunAnalysis')

class CheckStatus(Resource):
    def get(self, task_id):
        task = celery.AsyncResult(task_id)
        return {
            'id': task.id,
            'status': task.status,
            'result': task.result
        }, 200
api.add_resource(CheckStatus, '/check-status/<task_id>')

