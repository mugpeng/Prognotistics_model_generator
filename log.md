# 灵感

- 即使形如NC：Machine learning-based integration develops an immune-derived lncRNA signature for improving outcomes in colorectal cancer
  - 也就是基因筛选花了点功夫
  - 模型构建花了点功夫
  - 再就是集成学习、深度学习算法用来建模



# 230209





# 230208

- add filter input gene function(control the genes size below 50)
  - if length > 50, then add p value


- R parallels dont output when repeats dont finish

- timeout methods



- Run parallel in each model training steps
  - step.cox is too slow



# 230207

- run simple mode



# 230206

- model generator

- use gene symbol, mRNA(easy to find validation)
- diff mthods



## 数据

### 输入数据

tcga 的各个癌种，Msigdb 上的各种通路集。



### 特征数据

Msigdb N 种基因集。



## 方法

阈值设定。



### 特征筛选

Msigdb 上的各种基因集。



- 肿瘤正常差异基因
- ~~转移非转移差异基因~~
- cox 单因素生成相关基因
- lasso 筛选结果



- 与筛选基因集相关基因
- WGCNA 相关基因



### 模型选择

单一模型：

- lasso
- cox
- 随机森林
- GBM
- XGB
- NBC
- DT



## 其他

- 某癌种所有高ROC 结果的交集基因