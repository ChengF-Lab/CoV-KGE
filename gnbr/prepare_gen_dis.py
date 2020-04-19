import pandas as pd
from pandas import Series,DataFrame

pd.set_option("precision",10)
def  part1_pat2_concat():
    parti_gen_dis = pd.read_csv("original resource/part-i-gene-disease-path-theme-distributions.txt", delimiter="\t")
    print("done")
    parti_gen_dis.rename(columns={'path':'Dependence_path'},inplace=True)

    print(parti_gen_dis)
    # print(parti_gen_dis.describe())

    gen_dis_name = ['PubMed_ID', 'Sentence_number', 'Entity1', 'Loc1','Entity2','Loc2','Entity_Raw_str1','Entity_Raw_str2','DB_ID1','DB_ID2','Entity1_type','Entity2_type','Dependence_path','Sentence']
    print("name")
    partii_gen_dis =pd.read_csv("original resource/part-ii-dependency-paths-gene-disease-sorted-with-themes.txt", delimiter="\t", names = gen_dis_name, dtype={'DEVICE_ADDRESS': 'str'})
    print(partii_gen_dis)
    partii_gen_dis = partii_gen_dis[['Entity1','Entity2','DB_ID1','DB_ID2','Dependence_path']]

    print("过滤")
    #将依赖路劲改成小写：
    partii_gen_dis['Dependence_path'] = partii_gen_dis['Dependence_path'].map(lambda  x: str(x).lower())
    partii_gen_dis['Dependence_path'] = partii_gen_dis['Dependence_path'].map(lambda  x: str(x).lower())
    print('小写完成')

    #将两个表合并
    he_gen_dis = pd.merge(partii_gen_dis, parti_gen_dis, how='left', on='Dependence_path')
    print("合并完成")

    # 过滤合并后的表格U	U.ind	Ud	Ud.ind	D	D.ind	J	J.ind	Te	Te.ind	Y	Y.ind	G	G.ind	Md	Md.ind	X	X.ind	L	L.ind
    T_require = he_gen_dis['U'].map(lambda x: pd.notnull(x))
    C_require = he_gen_dis['Ud'].map(lambda x: pd.notnull(x))
    Sa_require = he_gen_dis['D'].map(lambda x: pd.notnull(x))
    Pr_require =he_gen_dis['J'].map(lambda x: pd.notnull(x))
    Pa_require =he_gen_dis['Te'].map(lambda x: pd.notnull(x))
    J_require =he_gen_dis['Y'].map(lambda x: pd.notnull(x))
    Mp_require =he_gen_dis['G'].map(lambda x: pd.notnull(x))
    J_require1=he_gen_dis['Md'].map(lambda x: pd.notnull(x))
    Mp_require1 =he_gen_dis['X'].map(lambda x: pd.notnull(x))
    Mp_require2 =he_gen_dis['L'].map(lambda x: pd.notnull(x))
    some = he_gen_dis[T_require |C_require| Sa_require | Pr_require | Pa_require|J_require|Mp_require|J_require1|Mp_require1|Mp_require2]


    #删除DBID为空的对
    Db1_require = some['DB_ID1'].map(lambda x: x!='null')
    Db2_require = some['DB_ID2'].map(lambda x: x!='null')

    triples_gen_dis = some[Db1_require & Db2_require]
    print("过滤合并完成")
    triples_gen_dis.to_csv("triples/triples_gen_dis_themes.tsv", sep="\t", header=True, index=False)

    #提取出实体,并修改实体DB名称
    triples_gen_dis ["DB_ID1"] = triples_gen_dis ["DB_ID1"] .apply(lambda x :"G:"+str(x))
    triples_gen_dis ["DB_ID2"] = triples_gen_dis ["DB_ID2"] .apply(lambda x :"D:"+str(x))

    gen_dis_gene = triples_gen_dis[["Entity1", 'DB_ID1']].drop_duplicates(['DB_ID1'])
    gen_dis_disease = triples_gen_dis[["Entity2", 'DB_ID2']].drop_duplicates(['DB_ID2'])

    gen_dis_gene.to_csv("original entitys/entity_gd_gene.tsv",sep='\t',header=True,index=False)
    print("to_csv chemical")
    gen_dis_disease.to_csv("original entitys/entity_gd_disease.tsv",sep='\t',header=True,index=False)
    print("entity done")

def relation_normalization():
    # #提取三元组
    triple_gen_dis= pd.read_csv("triples/triples_gen_dis_themes.tsv",delimiter='\t')
    print(triple_gen_dis.describe())

    print("最大值：",max(triple_gen_dis[['L']].values))
    #68514
    #将关系进行归一化：
    #U	U.ind	Ud	Ud.ind	D	D.ind	J	J.ind	Te	Te.ind	Y	Y.ind	G	G.ind	Md	Md.ind	X	X.ind	L	L.ind

    triple_gen_dis ["U"] = triple_gen_dis ["U"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["Ud"] = triple_gen_dis ["Ud"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["D"] = triple_gen_dis ["D"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["J"] = triple_gen_dis ["J"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["Te"] = triple_gen_dis ["Te"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["Y"] = triple_gen_dis ["Y"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["G"] = triple_gen_dis ["G"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["Md"] = triple_gen_dis ["Md"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["X"] = triple_gen_dis ["X"] .apply(lambda x : x/68514.0)
    triple_gen_dis ["L"] = triple_gen_dis ["L"] .apply(lambda x : x/68514.0)
    print("归一化完成")
    triples= []
    for index ,row  in triple_gen_dis.iterrows():
        if index %50000 ==0:
            print("读入",index/50000,"行")
        entity_1 = row["DB_ID1"]
        entity_2 = row["DB_ID2"]
        U = row['U']
        Ud = row['Ud']
        D = row['D']
        J = row['J']
        Te = row['Te']
        Y = row['Y']
        G = row['G']
        Md = row['Md']
        X = row['X']
        L = row['L']
        if U!=0.0:
            triples.append([entity_1,"U",entity_2,U])
        if Ud != 0.0:
            triples.append([entity_1, "Ud", entity_2, Ud])
        if D != 0.0:
            triples.append([entity_1, "D", entity_2, D])
        if J != 0.0:
            triples.append([entity_1, "J", entity_2, J])
        if Te != 0.0:
            triples.append([entity_1, "Te", entity_2, Te])
        if Y != 0.0:
            triples.append([entity_1, "Y", entity_2, Y])
        if G != 0.0:
            triples.append([entity_1, "G", entity_2, G])
        if Md != 0.0:
            triples.append([entity_1, "Md", entity_2, Md])
        if X != 0.0:
            triples.append([entity_1, "X", entity_2, X])
        if L != 0.0:
            triples.append([entity_1, "L", entity_2, L])
    print("read done")
    gen_dis_triples = DataFrame(triples,columns=["gene", "rel", "dis", "score"])

    gen_dis_triples.drop_duplicates(["gene", "rel", "dis", "score"],inplace=True)

    final_gen_dis_triples = gen_dis_triples.groupby(["gene", "rel", "dis"]).apply(lambda x: x.sort_values('score', ascending=False)).groupby(
        ["gene", "rel", "dis"]).first().reset_index()

    final_gen_dis_triples.to_csv("triples_gen_rel_dis.tsv",sep='\t',header=False,index=False)
    print(gen_dis_triples)


if __name__ == '__main__':
    part1_pat2_concat()
    relation_normalization()