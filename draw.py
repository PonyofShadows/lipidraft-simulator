import pymongo
import matplotlib
matplotlib.use('Agg') #非交互模式
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 150  # DPI默认为100
# 连接 MongoDB 数据库
client = pymongo.MongoClient("mongodb://localhost:27017/")

# 指定要遍历的数据库名
db_name = "test_lb_plus"

# 获取数据库对象
db = client[db_name]

# 遍历所有集合
for i in range(64):
    # 构建集合名
    collection_name = "exp_00_"+str(i)
    
    # 获取集合对象
    collection = db[collection_name]
    
    # 遍历集合中的所有文档
    dociter = collection.find()
    doc = dociter.next()
    # 处理每个文档数据    
    # 获取参数
    para = {}
    para["temp"] = doc["parameters"]["temperature"]
    para["ca"] = doc["parameters"]["molar fraction of A"]
    para["oc"] = doc["parameters"]['ccupancy fraction of C']
    para['la'] = doc['parameters']['tail chain length of A']
    para['lb'] = doc['parameters']['tail chain length of B']
    # 获取数据
    data = {
            'begin':{
                'phol':{},"chol":{}
                },
            'end':{
                'phol':{},"chol":{}
                }
            }
    for step in ['begin', 'end']:
        doc = dociter.next()
        para['energy'+step] = doc[step]['energy']
        data[step]['phol'] = doc[step]['phospholipids']
        data[step]['chol'] = doc[step]['cholesterol']
    #迭代次数
    doc = dociter.next()
    para['iteration'] = doc['iteration']
    #画图
    for step in ['begin','end']:
        #区分开A和B
        pldot={
                'A':{'x':[],'y':[]},
                'B':{'x':[],'y':[]}
            }
        dotnameslist = data[step]['phol']['phospholipid type']
        dotxs = iter(data[step]['phol']['x'])
        dotys = iter(data[step]['phol']['y'])
        for dotname in dotnameslist:
            if dotname == 'A':
                pldot['A']['x'].append(dotxs.__next__())
                pldot['A']['y'].append(dotys.__next__())
            elif dotname == 'B':
                pldot['B']['x'].append(dotxs.__next__())
                pldot['B']['y'].append(dotys.__next__())
        plt.scatter(
                pldot['A']['x'],
                pldot['A']['y'],
                color='red',
                s=8
                      )
        plt.scatter(
                pldot['B']['x'],
                pldot['B']['y'],
                color='blue',
                s=8
                    )
        plt.scatter(
                data[step]['chol']['x'],
                data[step]['chol']['y'],
                color='gray',
                s=2
                )
        plt.title(
                ('Lipid raft lattice of '+
                step+'\n'+
                 'ca:'+'{:.2f}'+'  '+
                 'lb:'+'{:.2f}'+'  '+
                 'oc:'+'{:.2f}'+'  ').format(para['ca'],para['lb'],para['oc'])
                )
        plt.xlabel('x (nm)')
        plt.ylabel('y (nm)')
        plt.savefig('./fig_lb_plus/'+'fig_'+str(i)+'_'+step)
        plt.clf()
