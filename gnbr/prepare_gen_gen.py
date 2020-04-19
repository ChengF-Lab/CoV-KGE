import pandas as pd
from pandas import Series,DataFrame

pd.set_option("precision",10)
def part1_pat2_concat():
    parti_gen_gen = pd.read_csv("original resource/part-i-gene-gene-path-theme-distributions.txt", delimiter="\t")
    print("done")
    print(parti_gen_gen)
    parti_gen_gen.rename(columns={'path':'Dependence_path'},inplace=True)

    gen_gen_name = ['PubMed_ID', 'Sentence_number', 'Entity1', 'Loc1','Entity2','Loc2','Entity_Raw_str1','Entity_Raw_str2','DB_ID1','DB_ID2','Entity1_type','Entity2_type','Dependence_path','Sentence']
    print("name")
    partii_gen_gen =pd.read_csv("original resource/part-ii-dependency-paths-gene-gene-sorted-with-themes.txt", delimiter="\t", names = gen_gen_name, dtype={'DEVICE_ADDRESS': 'str'})
    print("done")
    new_partii_gen_gen = partii_gen_gen[['Entity1','Entity2','DB_ID1','DB_ID2','Dependence_path']]
    print(new_partii_gen_gen)
    print("过滤")
    #将依赖路劲改成小写：
    new_partii_gen_gen['Dependence_path'] = new_partii_gen_gen['Dependence_path'].map(lambda  x: str(x).lower())
    parti_gen_gen['Dependence_path'] = parti_gen_gen['Dependence_path'].map(lambda  x: str(x).lower())
    print('小写完成')

    #将两个表合并
    he_gen_gen = pd.merge(new_partii_gen_gen, parti_gen_gen, how='left', on='Dependence_path')
    print("合并完成")

    #过滤合并后的表格

    T_require = he_gen_gen['B'].map(lambda x: pd.notnull(x))
    C_require = he_gen_gen['W'].map(lambda x: pd.notnull(x))
    Sa_require = he_gen_gen['V+'].map(lambda x: pd.notnull(x))
    Pr_require =he_gen_gen['E+'].map(lambda x: pd.notnull(x))
    Pa_require =he_gen_gen['E'].map(lambda x: pd.notnull(x))
    J_require =he_gen_gen['I'].map(lambda x: pd.notnull(x))
    Mp_require =he_gen_gen['H'].map(lambda x: pd.notnull(x))
    J_require1=he_gen_gen['Rg'].map(lambda x: pd.notnull(x))
    Mp_require1 =he_gen_gen['Q'].map(lambda x: pd.notnull(x))
    some = he_gen_gen[T_require |C_require| Sa_require | Pr_require | Pa_require|J_require|Mp_require|J_require1|Mp_require1]
    print("过滤为空的主题")
    #
    #删除DBID为空的对
    Db1_require = some['DB_ID1'].map(lambda x: x!='null')
    Db2_require = some['DB_ID2'].map(lambda x: x!='null')

    triples_gen_gen = some[Db1_require & Db2_require]
    #
    print("过滤合并完成")

    #
    #提取出实体,并修改实体DB名称
    triples_gen_gen ["DB_ID1"] = triples_gen_gen ["DB_ID1"] .apply(lambda x :"G:"+str(x))
    triples_gen_gen ["DB_ID2"] = triples_gen_gen ["DB_ID2"] .apply(lambda x :"G:"+str(x))

    new_triples_gen_gen = triples_gen_gen[["Entity1","DB_ID1","Entity2","DB_ID2",'B','B.ind','W','W.ind','V+','V+.ind','E+','E+.ind','E','E.ind','I','I.ind','H','H.ind','Rg','Rg.ind','Q','Q.ind']]
    new_triples_gen_gen.to_csv("triples/triples_gen_gen_themes.tsv",sep="\t",header=True)
    print(new_triples_gen_gen)

    gen_gen_gene1= triples_gen_gen[["Entity1",'DB_ID1']].drop_duplicates(['DB_ID1'])
    gen_gen_gene2 = triples_gen_gen[["Entity2",'DB_ID2']].drop_duplicates(['DB_ID2'])
    print("gene1:")
    print(gen_gen_gene1)
    print("gene2：")
    print(gen_gen_gene2)
    gen_gen_gene1.to_csv("original entitys/entity_gg_gene1.tsv",sep="\t",header=True,index=False)
    gen_gen_gene2.to_csv("original entitys/entity_gg_gene2.tsv",sep="\t",header=True,index=False)
#提取三元组

# relation2id = pd.read_csv("Relation2Id.tsv",delimiter='\t')
#
def relation_normalization():
    triple_gen_gen= pd.read_csv("triples/triples_gen_gen_themes.tsv",delimiter='\t')
    print(triple_gen_gen.describe())
    max = max(triple_gen_gen[['Q']].values)
    print("max values:",max)
    #515159
    #将关系进行归一化：
    #B	B.ind	W	W.ind	V+	V+.ind	E+	E+.ind	E	E.ind	I	I.ind	H	H.ind	Rg	Rg.ind	Q	Q.ind

    triple_gen_gen ["B"] = triple_gen_gen ["B"] .apply(lambda x : x/max)
    triple_gen_gen ["W"] = triple_gen_gen ["W"] .apply(lambda x : x/max)
    triple_gen_gen ["V+"] = triple_gen_gen ["V+"] .apply(lambda x : x/max)
    triple_gen_gen ["E+"] = triple_gen_gen ["E+"] .apply(lambda x : x/max)
    triple_gen_gen ["E"] = triple_gen_gen ["E"] .apply(lambda x : x/max)
    triple_gen_gen ["I"] = triple_gen_gen ["I"] .apply(lambda x : x/max)
    triple_gen_gen ["H"] = triple_gen_gen ["H"] .apply(lambda x : x/max)
    triple_gen_gen ["Rg"] = triple_gen_gen ["Rg"] .apply(lambda x : x/max)
    triple_gen_gen ["Q"] = triple_gen_gen ["Q"] .apply(lambda x : x/max)

    print("归一化完成")
    triples= []
    for index ,row  in triple_gen_gen.iterrows():
        if index %50000 ==0:
            print("读入",index/50000,"行")
        entity_1 = row["DB_ID1"]
        entity_2 = row["DB_ID2"]
        B = row['B']
        W = row['W']
        V1 = row['V+']
        E1 = row['E+']
        E = row['E']
        I = row['I']
        H = row['H']
        Rg = row['Rg']
        Q = row['Q']

        if B!=0.0:
            triples.append([entity_1,"B",entity_2,B])
        if W != 0.0:
            triples.append([entity_1, "W", entity_2, W])
        if V1 != 0.0:
            triples.append([entity_1, "V+", entity_2, V1])
        if E1 != 0.0:
            triples.append([entity_1, "E+", entity_2, E1])
        if E != 0.0:
            triples.append([entity_1, "E", entity_2, E])
        if I != 0.0:
            triples.append([entity_1, "I", entity_2, I])
        if H != 0.0:
            triples.append([entity_1, "H", entity_2, H])
        if Rg != 0.0:
            triples.append([entity_1, "Rg", entity_2, Rg])
        if Q != 0.0:
            triples.append([entity_1, "Q", entity_2, Q])

    print("read done")
    gen_gen_triples = DataFrame(triples,columns=["gene", "rel", "gen", "score"])

    gen_gen_triples.drop_duplicates(["gene", "rel", "gen", "score"],inplace=True)
    final_gen_gen_triples= gen_gen_triples.sort_values('score', ascending=False).groupby(["ent1", "rel", "ent2"]).first().reset_index()
    final_gen_gen_triples.to_csv("triples_gen_rel_gen.tsv",sep='\t',header=False,index=False)
    print(gen_gen_triples)

if __name__ == '__main__':
    part1_pat2_concat()
    relation_normalization()
