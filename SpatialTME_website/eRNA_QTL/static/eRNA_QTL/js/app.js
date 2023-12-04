"use strict";
var app = angular.module('eRNA_QTL', ['ui.bootstrap', 'ngRoute','ui.bootstrap-slider']);
// ,'chart.js'
    app.config(function ($routeProvider) {
        $routeProvider
            .when("/", {
                templateUrl: "/static/eRNA_QTL/pages/home.html",
                controller: "HomeController",
            })
            .when("/browse", {
                templateUrl: "/static/eRNA_QTL/pages/browse.html",
                controller: "BrowseController",
            })
            .when("/search", {
                templateUrl: "/static/eRNA_QTL/pages/search.html",
                controller: "SearchController",
            })
            .when("/browse/:cancer_type", {
                templateUrl: "/static/eRNA_QTL/pages/cancer_type.html",
                controller: "CancerTypeController",
            })
            .when("/search/:platform/:dataset_name/:sample_id", {
                templateUrl: "/static/eRNA_QTL/pages/analysis.html",
                controller: "AnalysisController",
            })
            .when("/analysis/:module", {
                templateUrl: "/static/eRNA_QTL/pages/module_analysis.html",
                controller: "ModuleAnalysisController",
            })
            .when("/upload", {
                templateUrl: "/static/eRNA_QTL/pages/upload.html",
                controller: "UploadController",
            })
            .when("/tutorial", {
                templateUrl: "/static/eRNA_QTL/pages/tutorial.html",
                controller: "TutorialController",
            })
            .when("/test", {
                templateUrl: "/static/eRNA_QTL/pages/test.html",
                controller: "TestController",
            })
            .otherwise({
                redirectTo: "/404.html",
            });
    })
    app.config(function ($interpolateProvider) {
        $interpolateProvider.startSymbol('{$');
        $interpolateProvider.endSymbol('$}');
    })
    app.service('eRNA_QTLService',function () {
        this.getAPIBaseUrl = function () {
             return "";
            //return "/eRNA_QTL/";
        }
    })
    app.filter('unique', function () {
        return function (collection, keyname) {
          var output = [],
            keys = [];
          angular.forEach(collection, function (item) {
            var key = item[keyname];
            if (keys.indexOf(key) === -1) {
              keys.push(key);
              output.push(item);
            }
          });
          return output;
        };
      });
    app.filter('trustAsResourceUrl', ['$sce', function($sce) {
        return function(val) {
            return $sce.trustAsResourceUrl(val);
        };
    }])

angular.module('eRNA_QTL')
    .controller('BrowseController', BrowseController);
function BrowseController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    console.log("BrowseController loaded");
    $window.scrollTo(0, 0);

    var base_url = eRNA_QTLService.getAPIBaseUrl()

    $scope.records = [
        { organ:'Breast', pic_src: "/static/eRNA_QTL/image/breast.png", 
         cancerlist_num: ['BRCA (59)'],
         cancerlist:['BRCA']},
        { organ:'Colorectum', pic_src: "/static/eRNA_QTL/image/colon.png", 
        cancerlist_num: ['CRC (14)'],
        cancerlist:['CRC']},
        { organ:'Liver', pic_src: "/static/eRNA_QTL/image/liver.png", 
        cancerlist_num: ['HCC (27)'],
        cancerlist: ['HCC']},
        { organ:'Kidney', pic_src: "/static/eRNA_QTL/image/kidney.png", 
        cancerlist_num: ['RCC (24)'],
        cancerlist: ['RCC']},
         { organ:'Pancreas', pic_src: "/static/eRNA_QTL/image/pancreas.png", 
         cancerlist_num: ['PDAC (4)'],
         cancerlist: ['PDAC']},
         { organ:'Prostate', pic_src: "/static/eRNA_QTL/image/prostate.png", 
         cancerlist_num: ['PRAD (24)'],
         cancerlist: ['PRAD']},
         { organ:'Stomach', pic_src: "/static/eRNA_QTL/image/stomach.png", 
         cancerlist_num: ['GIST (2)'],
         cancerlist: ['GIST']},
         { organ:'Lung', pic_src: "/static/eRNA_QTL/image/lung.png", 
         cancerlist_num: ['LUSC (1)','LUAD (3)'],
         cancerlist: ['LUSC','LUAD']},
         { organ:'Head', pic_src: "/static/eRNA_QTL/image/head.png", 
         cancerlist_num: ['HNSC (4)'],
         cancerlist: ['HNSC']},
         { organ:'Skin', pic_src: "/static/eRNA_QTL/image/skin.png", 
         cancerlist_num: ['SCC (8)','SKCM (2)'],
         cancerlist: ['SCC','SKCM']},
         { organ:'Pelvic cavity', pic_src: "/static/eRNA_QTL/image/pelvic.png", 
         cancerlist_num: ['OV (37)',
                      'CSCC (1)',
                      'UCEC (1)'],
         cancerlist: ['OV',
                      'CSCC',
                      'UCEC']},
        { organ:'Brain', pic_src: "/static/eRNA_QTL/image/brain.png", 
        cancerlist_num: ['GBM (32)',
                      'EPN (14)',
                      'MB (3)',
                      'PCNSL (4)'],
        cancerlist: [ 'GBM',
                      'EPN',
                      'MB',
                      'PCNSL']
        }
        
      ];
}
angular.module('eRNA_QTL')
    .controller('AnalysisController', AnalysisController);

function AnalysisController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    $window.scrollTo(0, 0);
    console.log("AnalysisController loaded");
    $('.selectpicker').selectpicker('refresh');
    var base_url = eRNA_QTLService.getAPIBaseUrl()
    var sample_id = $routeParams.sample_id
    var dataset_name = $routeParams.dataset_name
    var platform = $routeParams.platform

    console.log($routeParams.sample_id);
    // var analysis_URL="http://139.196.214.89:3838/SpatialTME_shiny/#!/?dataset="+dataset_name+"&sample="+sample_id
    var analysis_URL="http://139.196.214.89:3838/SpatialTME_shiny/#!/?platform="+platform+"&dataset="+dataset_name+"&sample="+sample_id

    $scope.analysis_URL=analysis_URL

    $http({
        url: base_url + '/api/samples',
        method: 'GET',
        params: {'organ_flask':'',
                'cancer_type_flask':'',
                'accessible_ID_flask':'',
                'platform':'',
                'sample_id_flask':$routeParams.sample_id}, 
        }).then(
        function (response) {
            var Sample_list = response.data.Sample_list;
            $scope.sample = Sample_list
            console.log(Sample_list);
        }
    )
    $scope.getSampleImageURL = function(sample) {
        var imageURL = '/static/eRNA_QTL/HE/' + sample.sample_id + '_tissue_lowres_image.png';
        return imageURL;
      };

}

angular.module('eRNA_QTL')
    .controller('CancerTypeController', CancerTypeController);

function CancerTypeController($scope,$http,$window,$routeParams,eRNA_QTLService,$location) {
    console.log("CancerTypeController loaded");
    $window.scrollTo(0, 0);
    $('.selectpicker').selectpicker('refresh');
    var base_url = eRNA_QTLService.getAPIBaseUrl()
    var cancer_type = $routeParams.cancer_type
    $scope.cancer_type = cancer_type
    $http({
        url: base_url + '/api/samples',
        method: 'GET',
        params:  {'cancer_type_flask':cancer_type,'organ_flask':"",'accessible_ID_flask':'','sample_id_flask':'','platform':''},
    }).then(
        function (response) {
            var Sample_list = response.data.Sample_list;
            $scope.Sample_list = Sample_list
            $(document).ready(function () {
                $('#tablesort_cancertype').DataTable(
                    {
                        scrollX: true
                    });
            })
        }
    )
}

angular.module('eRNA_QTL')
    .controller('HomeController', HomeController);

function HomeController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    console.log("HomeController loaded");
    $window.scrollTo(0, 0);

    var base_url = eRNA_QTLService.getAPIBaseUrl()
}

angular.module('eRNA_QTL')
    .controller('ModuleAnalysisController', ModuleAnalysisController);

function ModuleAnalysisController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    console.log("ModuleAnalysisController loaded");
    $window.scrollTo(0, 0);
    $('.selectpicker').selectpicker('refresh');
    var base_url = eRNA_QTLService.getAPIBaseUrl()
    var module = $routeParams.module
    $scope.module = module
    console.log(module);
    $scope.data=[]
    $scope.cancerFilter = ''; // 用户选择的 cancer type
    $scope.sampleFilter = ''; // 用户选择的 sample
    $scope.analysis_result=0
    $scope.sample_info=0
    $scope.analysis_URL = '';


    $http({
        url: base_url + '/api/samples',
        method: 'GET',
        params: {'organ_flask':'',
                'cancer_type_flask':'',
                'accessible_ID_flask':'',
                'sample_id_flask':'',
                'platform':''}, 
        }).then(
        function (response) {
            var data = response.data.Sample_list;
            // 如果模块是 "Sub-spot GEP"，则只保留特定 Visium 的样本
            if ($scope.module === 'Sub-spot_GEP') {
                data = data.filter(sample => sample.platform === 'Visium');
            }

            $scope.data = data;

        }
    )

      $scope.clear = function() {
        // 清空所有条件选择框
        $scope.cancerFilter = '';
        $scope.sampleFilter = '';
        $scope.analysis_result=0;
        $scope.sample_info=0
      };

      $scope.analysis=function(){
        $scope.analysis_result=0
        $scope.sample_info=0
        var selectedSample = $scope.data.find(function(sample) {
            return sample.sample_id === $scope.sampleFilter.sample_id;
          });

        if (selectedSample) {
            var datasetName = selectedSample.dataset_name;
            // 在这里处理获取到的 datasetName
        }
        $scope.sample = selectedSample
        var analysis_URL="http://139.196.214.89:3838/SpatialTME_shiny/#!/?platform="+selectedSample.platform+"&dataset="+selectedSample.dataset_name+"&sample="+selectedSample.sample_id+"&selectedTab="+module

        console.log(analysis_URL);
        $scope.analysis_URL=analysis_URL
        $scope.sample_info=1
        $scope.analysis_result=1
    }

    $scope.getSampleImageURL = function(sample) {
        var imageURL = '/static/eRNA_QTL/HE/' + sample.sample_id + '_tissue_lowres_image.png';
        return imageURL;
      };
      
}

angular.module('eRNA_QTL')
    .controller('SearchController', SearchController);

function SearchController($scope,$http,$window,$routeParams,eRNA_QTLService,$location) {
    console.log("SearchController loaded");
    $window.scrollTo(0, 0);
    $('.selectpicker').selectpicker('refresh');
    var base_url = eRNA_QTLService.getAPIBaseUrl()
    
    $http({
        url: base_url + '/api/samples',
        method: 'GET',
        params:  {'cancer_type_flask':"",'organ_flask':"",'accessible_ID_flask':'','sample_id_flask':'','platform':''},
    }).then(
        function (response) {
            var Init_Sample_list = response.data.Sample_list;
            $scope.Init_Sample_list = Init_Sample_list
        }
    )
    $scope.cancerFilter="";
    $scope.organFilter="";
    $scope.projectFilter="";
    $scope.platformFilter = ""; // 初始化platformFilter

    $scope.search=function () {
        $scope.table_show=0
        $scope.loading=1

        $http({
            url: base_url + '/api/samples',
            method: 'GET',
            params:  {'organ_flask':$scope.organFilter,
                      'cancer_type_flask':$scope.cancerFilter,
                      'accessible_ID_flask':$scope.projectFilter,
                      'sample_id_flask':'',
                      'platform':$scope.platformFilter},
        }).then(
            function (response) {
                // init table
                $('a[data-toggle="tab"]').on( 'shown.bs.tab', function (e) {
                    $.fn.dataTable.tables( {visible: true, api: true} ).columns.adjust();
                } );
                $scope.load = 0
                $scope.Sample_list = response.data.Sample_list;
                $scope.tabletag=1
                $scope.table_show=1
                $scope.loading=0
                $('#tablesort').DataTable().destroy();// destroy exist datatable
                $(document).ready(function () {
                    $('#tablesort').DataTable(
                        {
                            destroy: true,//init datble
                            scrollX: true
                        });
                })
            }
        )
    }

    $scope.openInNewWindow = function(sample) {
        var url = '#!/search/' + sample.platform + '/' + sample.dataset_name + '/' + sample.sample_id;
        window.open(url, '_blank');
    };


    $scope.clear = function() {
        // 清空所有条件选择框
        $scope.organFilter = '';
        $scope.projectFilter = '';
        $scope.cancerFilter = '';
        $scope.platformFilter = ''
      };

}
angular.module('eRNA_QTL')
    .controller('TutorialController', TutorialController);

function TutorialController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    console.log("TutorialController loaded");
    $window.scrollTo(0, 0);

    var base_url = eRNA_QTLService.getAPIBaseUrl()
}

angular.module('eRNA_QTL')
    .controller('UploadController', UploadController);

function UploadController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    console.log("UploadController loaded");
    $window.scrollTo(0, 0);
    $scope.uploadSuccess = false;
    $scope.progressVisible = false;
    $scope.showRunButton = false; // 初始状态下隐藏Run按钮

    $scope.openModal = function openModal(modalId) {
        var modal = document.getElementById(modalId);
        modal.style.display = "block";
      }
    
    // JavaScript函数，关闭指定ID的弹窗
    $scope.closeModal = function closeModal(modalId) {
        var modal = document.getElementById(modalId);
        modal.style.display = "none";
    }

    // 为 ST 文件上传添加事件监听器
    $scope.stFileChanged = function(element) {
        $scope.$apply(function() {
            var file = element.files[0];
            $scope.stFileName = file ? file.name : '';
        });
    };
    // 为 scRNA 文件上传添加事件监听器
    $scope.scrnaMatFileChanged = function(element) {
        $scope.$apply(function() {
            var file = element.files[0];
            $scope.scrnaMatFileName = file ? file.name : '';
        });
    };
    $scope.scrnaMetadataFileChanged = function(element) {
        $scope.$apply(function() {
            var file = element.files[0];
            $scope.scrnaMetadataFileName = file ? file.name : '';
        });
    };
    
    document.getElementById('clear').addEventListener('click', function() {
        // 清空文件上传输入
        document.getElementById('st-file-upload').value = '';
        document.getElementById('scrna-mat-file-upload').value = '';
        document.getElementById('scrna-metadata-file-upload').value = '';
            // 清空文件名显示
        $scope.$apply(function() {
            $scope.stFileName = '';
            $scope.scrnaMatFileName = '';
            $scope.scrnaMetadataFileName = '';
        });
        // 清空匹配单元格类型的输入
        document.getElementById('malignant-cluster').value = '';
        document.getElementById('tissue-cluster').value = '';
        document.getElementById('stromal-cluster').value = '';
        document.getElementById('email').value = '';
      });
      
    // upload/
    document.getElementById('submit').addEventListener('click', function(e) {
        e.preventDefault();

        var stFile = document.getElementById('st-file-upload').files[0];
        var scrnaMatFile = document.getElementById('scrna-mat-file-upload').files[0];
        var scrnaMetadataFile = document.getElementById('scrna-metadata-file-upload').files[0];
        var tumorCell = document.getElementById('malignant-cluster').value;
        var epithelialCell = document.getElementById('tissue-cluster').value;
        var stromalCell = document.getElementById('stromal-cluster').value;
        var email = document.getElementById('email').value;
        var emailRegex = /^[a-zA-Z0-9._-]+@[a-zA-Z0-9.-]+\.(edu|ac|org|edu\.[a-zA-Z]{2,3}|org\.[a-zA-Z]{2,3}|ac\.[a-zA-Z]{2,3})$/;

        if (!stFile || !scrnaMatFile|| !scrnaMetadataFile || !tumorCell || !epithelialCell || !stromalCell|| !email) {
            alert('Please fill out all fields.');
            return;
        }
        
        if (!validateFileTypes(stFile, scrnaMatFile, scrnaMetadataFile)) {
            return;  // 如果验证失败，停止执行
        }

        if (!emailRegex.test(email)) {
            alert('Please enter a valid email address.');
            return; // 如果电子邮件格式不正确，则停止执行函数
        }

        // 检查文件大小
        var maxFileSize = 300 * 1024 * 1024; // 300MB
        if ((stFile && stFile.size > maxFileSize) || 
            (scrnaMatFile && scrnaMatFile.size > maxFileSize) || 
            (scrnaMetadataFile && scrnaMetadataFile.size > maxFileSize)) {
            alert('File size should not exceed 300MB.');
            return;
        }

        var formData = new FormData();
        formData.append('stFile', stFile);
        formData.append('scrnaMatFile', scrnaMatFile);
        formData.append('scrnaMetadataFile', scrnaMetadataFile);
        formData.append('tumorCell', tumorCell);
        formData.append('epithelialCell', epithelialCell);
        formData.append('stromalCell', stromalCell);
        formData.append('email', email);
        console.log(formData)

        var xhr = new XMLHttpRequest();
        
        xhr.upload.addEventListener("progress", function (evt) {
            if (evt.lengthComputable) {
                var percentComplete = Math.round((evt.loaded / evt.total) * 100);
              // 更新进度条的值
              document.getElementById('uploadProgress').value = percentComplete;
              $scope.progressVisible = true;
              $scope.$apply(); // 触发 $scope 的更新，以便反映在视图上
            }
        }, false);

        xhr.addEventListener("load", function () {
            if (xhr.status == 200) {
              console.log("Upload complete!");
              $scope.$apply(function () {
                $scope.uploadSuccess = true;
                $scope.progressVisible = false;
                $scope.showRunButton = true; // 显示Run按钮
              });
            }else{
                console.error("Upload failed!");
                // 显示从后端返回的错误消息
                var errorMsg = JSON.parse(xhr.responseText).message;
                alert(errorMsg);
                $scope.progressVisible = false;
            }
        }, false);  

        xhr.open("POST", '/FileUpload', true);
        xhr.send(formData);
    });
    function validateFileTypes(stFile, scrnaMatFile, scrnaMetadataFile) {
        if (!stFile || !scrnaMatFile || !scrnaMetadataFile) {
            alert('Please upload all required files.');
            return false;
        }
    
        if (!/\.(tar\.gz|tar|zip)$/.test(stFile.name)) {
            alert('ST file must be a .tar.gz, .tar, or .zip file');
            return false;
        }
    
        if (!/\.(rds|csv)(\.gz)?$/.test(scrnaMatFile.name)) {
            alert('scRNA matrix file must be a .rds, .csv, .rds.gz, or .csv.gz file');
            return false;
        }
        
        if (!/\.csv$/.test(scrnaMetadataFile.name)) {
            alert('scRNA metadata file must be a .csv file');
            return false;
        }
    
        return true;
    }
    // Run/
    var runButton = document.getElementById('run');
    var isRunning = false;  // 增加一个标志变量来追踪运行状态

    if (runButton) {
        runButton.addEventListener('click', function(){
            if (isRunning) {
                alert('A process is already running, please wait for the current process to complete.');
                return; // 如果已经在运行，则不执行后续代码
            }
    
            // 更新运行状态
            isRunning = true;
            console.log("Run button clicked"); // 输出日志

            var stFile = document.getElementById('st-file-upload').files[0];
            var scrnaMatFile = document.getElementById('scrna-mat-file-upload').files[0];
            var scrnaMetadataFile = document.getElementById('scrna-metadata-file-upload').files[0];
            var tumorCell = document.getElementById('malignant-cluster').value;
            var epithelialCell = document.getElementById('tissue-cluster').value;
            var stromalCell = document.getElementById('stromal-cluster').value;
            var email = document.getElementById('email').value;

            if (!stFile || !scrnaMatFile || !scrnaMetadataFile || !tumorCell || !epithelialCell || !stromalCell|| !email) {
                alert('Please upload ST and scRNA data and fill out all fields.');
            } else {
                alert("The R script is running in the background, and the results will be sent to your email."); // 立即显示弹窗
                var formData = new FormData();
                formData.append('stFile', stFile);
                formData.append('scrnaMatFile', scrnaMatFile);
                formData.append('scrnaMetadataFile', scrnaMetadataFile);
                formData.append('tumorCell', tumorCell);
                formData.append('epithelialCell', epithelialCell);
                formData.append('stromalCell', stromalCell);
                formData.append('email', email);
    
                var xhr = new XMLHttpRequest();
                            
                xhr.addEventListener("load", function () {
                    if (xhr.status == 200) {
                        console.log("Run complete!");
                        isRunning = false;  // 重置运行状态
                    } else if (xhr.status === 404) {
                        alert("The corresponding upload folder was not found."); // 弹出窗口提示用户
                        isRunning = false;  // 重置运行状态
                    } else {
                        console.error("Run failed!");
                        console.error(xhr.responseText);
                        isRunning = false;  // 重置运行状态
                    }
                }, false);
        
                xhr.open("POST", '/RunAnalysis', true);
                xhr.send(formData);
            }
        });
    }else{
        console.error("Run button not found");
    }

}

angular.module('eRNA_QTL')
    .controller('TestController', TestController);

function TestController($scope,$http,$window,$routeParams,eRNA_QTLService) {
    console.log("TestController loaded");
    $window.scrollTo(0, 0);
    $scope.uploadSuccess = false;
    $scope.progressVisible = false;

    $scope.fileChanged = function(element) {
        $scope.$apply(function(scope) {
            $scope.file = element.files[0];
            console.log($scope.file);
        });
    };
    $scope.upload = function () {
        var formData = new FormData();
        formData.append('file',$scope.file);
        var xhr = new XMLHttpRequest();

        xhr.upload.addEventListener("progress", function (evt) {
            if (evt.lengthComputable) {
                var percentComplete = Math.round((evt.loaded / evt.total) * 100);
              // 更新进度条的值
              document.getElementById('uploadProgress').value = percentComplete;
              $scope.progressVisible = true;
              $scope.$apply(); // 触发 $scope 的更新，以便反映在视图上
            }
        }, false);

        xhr.addEventListener("load", function () {
            if (xhr.status == 200) {
              console.log("Upload complete!");
              $scope.$apply(function () {
                $scope.uploadSuccess = true;
                $scope.progressVisible = false;
              });
            } else {
              console.error("Upload failed!");
            }
        }, false);  
        xhr.open("POST", '/FileUpload', true);
        xhr.send(formData);
    };
}