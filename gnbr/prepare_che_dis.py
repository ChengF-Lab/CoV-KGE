import pandas as pd
from pandas import Series,DataFrame
pd.set_option("precision",10)
def part1_pat2_concat():
    parti_che_dis = pd.read_csv("original resource/part-i-chemical-disease-path-theme-distributions.txt",
                                delimiter="\t")
    print("done")
    parti_che_dis.rename(columns={'path': 'Dependence_path'}, inplace=True)

    print(parti_che_dis)
    # max
    che_dis_name = ['PubMed_ID', 'Sentence_number', 'Entity1', 'Loc1', 'Entity2', 'Loc2', 'Entity_Raw_str1',
                    'Entity_Raw_str2', 'DB_ID1', 'DB_ID2', 'Entity1_type', 'Entity2_type', 'Dependence_path',
                    'Sentence']
    partii_che_dis = pd.read_csv("original resource/part-ii-dependency-paths-chemical-disease-sorted-with-themes.txt",
                                 delimiter="\t", names=che_dis_name, dtype={'DEVICE_ADDRESS': 'str'})
    print("done")

    partii_che_dis = partii_che_dis[['Entity1', 'Entity2', 'DB_ID1', 'DB_ID2', 'Dependence_path']]

    #将依赖路劲改成小写：
    parti_che_dis['Dependence_path'] = parti_che_dis['Dependence_path'].map(lambda  x: x.lower())
    partii_che_dis['Dependence_path'] = partii_che_dis['Dependence_path'].map(lambda  x: x.lower())
    print('lower done')

    #merge two tables
    he_che_gen = pd.merge(partii_che_dis, parti_che_dis, how='left', on='Dependence_path')
    print("merge done")


    # filter null relation
    T_require = he_che_gen['T'].map(lambda x: pd.notnull(x))
    C_require = he_che_gen['C'].map(lambda x: pd.notnull(x))
    Sa_require = he_che_gen['Sa'].map(lambda x: pd.notnull(x))
    Pr_require =he_che_gen['Pr'].map(lambda x: pd.notnull(x))
    Pa_require =he_che_gen['Pa'].map(lambda x: pd.notnull(x))
    J_require =he_che_gen['J'].map(lambda x: pd.notnull(x))
    Mp_require =he_che_gen['Mp'].map(lambda x: pd.notnull(x))

    some = he_che_gen[T_require |C_require| Sa_require | Pr_require | Pa_require|J_require|Mp_require]
    print(some.describe())

    #delete DBID is null
    Db1_require = some['DB_ID1'].map(lambda x: x!='null')
    Db2_require = some['DB_ID2'].map(lambda x: x!='null')

    triples_che_dis = some[Db1_require & Db2_require]
    print("filter done")

    triples_che_dis.to_csv("triples/triples_che_dis_themes.tsv",sep='\t',header=True,index=False)

    # get entities and rename it C:chemical D:disease
    triples_che_dis ["DB_ID1"] = triples_che_dis ["DB_ID1"] .apply(lambda x :"C:"+str(x))
    triples_che_dis ["DB_ID2"] = triples_che_dis ["DB_ID2"] .apply(lambda x :"D:"+str(x))

    triples_che_dis.to_csv("triples/triples_che_dis_themes.tsv",sep="\t",header=True,index=False)
    print(triples_che_dis)
    #drop duplicate entities
    che_dis_chemical= triples_che_dis[["Entity1",'DB_ID1']].drop_duplicates(['DB_ID1'])
    che_dis_disease = triples_che_dis[["Entity2",'DB_ID2']].drop_duplicates(['DB_ID2'])

    print(che_dis_chemical)
    print(che_dis_disease)
    che_dis_chemical.to_csv("original entitys/entity_cd_chemical.tsv",sep='\t',header=True,index=False)
    print("to_csv chemical")
    che_dis_disease.to_csv("original entitys/entity_cd_disease.tsv",sep='\t',header=True,index=False)
    print("to_csv disease")


#形成三元组
def relation_normalization():

    triple_che_dis = pd.read_csv("triples/triples_che_dis_themes.tsv",delimiter='\t')
    print(triple_che_dis.describe())

    print("最大值：",max(triple_che_dis[['T']].values),max(triple_che_dis[['C']].values))
    #248178
    #将关系进行归一化：
    triple_che_dis ["T"] = triple_che_dis ["T"] .apply(lambda x : x/248178.0)
    triple_che_dis ["C"] = triple_che_dis ["C"] .apply(lambda x : x/248178.0)
    triple_che_dis ["Sa"] = triple_che_dis ["Sa"] .apply(lambda x : x/248178.0)
    triple_che_dis ["Pr"] = triple_che_dis ["Pr"] .apply(lambda x : x/248178.0)
    triple_che_dis ["Pa"] = triple_che_dis ["Pa"] .apply(lambda x : x/248178.0)
    triple_che_dis ["J"] = triple_che_dis ["J"] .apply(lambda x : x/248178.0)
    triple_che_dis ["Mp"] = triple_che_dis ["Mp"] .apply(lambda x : x/248178.0)
    print("归一化完成")
    #每一个不为0的主题都成为成为三元组
    triples= []
    def load():
        for index ,row  in triple_che_dis.iterrows():
            if index %50000 ==0:
                print("读入",index/50000,"行")
            entity_1 = row["DB_ID1"]
            entity_2 = row["DB_ID2"]
            T = row['T']
            C = row['C']
            Sa = row['Sa']
            Pr = row['Pr']
            Pa = row['Pa']
            J = row['J']
            Mp = row['Mp']
            if T!=0.0:
                triples.append([entity_1,"T",entity_2,T])
            if C != 0.0:
                triples.append([entity_1, "C", entity_2, C])
            if Sa != 0.0:
                triples.append([entity_1, "Sa", entity_2, Sa])
            if Pr != 0.0:
                triples.append([entity_1, "Pr", entity_2, Pr])
            if Pa != 0.0:
                triples.append([entity_1, "Pa", entity_2, Pa])
            if J != 0.0:
                triples.append([entity_1, "J", entity_2, J])
            if Mp != 0.0:
                triples.append([entity_1, "Mp", entity_2, Mp])
        print("read done")
        che_dis_triples = DataFrame(triples,columns=["che","rel","dis","score"])

        #去重：
        che_dis_triples.drop_duplicates(["che","rel","dis","score"],inplace=True)

        final_che_dis_triples= che_dis_triples.sort_values('score', ascending=False).groupby(["ent1", "rel", "ent2"]).first().reset_index()
        final_che_dis_triples.to_csv("triples_che_rel_dis.tsv",sep='\t',header=False,index=False)
        print(final_che_dis_triples)

if __name__ == '__main__':
    part1_pat2_concat()
    relation_normalization()