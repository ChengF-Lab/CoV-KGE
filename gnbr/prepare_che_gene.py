import pandas as pd
from pandas import Series,DataFrame
pd.set_option("precision",10)
def part1_pat2_concat():
    parti_che_gene = pd.read_csv("original resource/part-i-chemical-gene-path-theme-distributions.txt", delimiter="\t")
    print("done")
    parti_che_gene.rename(columns={'path':'Dependence_path'},inplace=True)

    print(parti_che_gene)
    #max 45222
    che_dis_name = ['PubMed_ID', 'Sentence_number', 'Entity1', 'Loc1','Entity2','Loc2','Entity_Raw_str1','Entity_Raw_str2','DB_ID1','DB_ID2','Entity1_type','Entity2_type','Dependence_path','Sentence']
    partii_che_gene =pd.read_csv("original resource/part-ii-dependency-paths-chemical-gene-sorted-with-themes.txt", delimiter="\t", names = che_dis_name, dtype={'DEVICE_ADDRESS': 'str'})
    print("done")
    partii_che_gene = partii_che_gene[['Entity1','Entity2','DB_ID1','DB_ID2','Dependence_path']]
    print(partii_che_gene)

    #将依赖路劲改成小写：
    parti_che_gene['Dependence_path'] = parti_che_gene['Dependence_path'].map(lambda  x: str(x).lower())
    partii_che_gene['Dependence_path'] = partii_che_gene['Dependence_path'].map(lambda  x: str(x).lower())
    print('小写完成')

    #将两个表合并
    he_che_gen = pd.merge(partii_che_gene, parti_che_gene, how='left', on='Dependence_path')
    print("合并完成")
    # print(he_che_gen.describe())
    # print(he_che_gen[1:10])

    #过滤合并后的表格

    T_require = he_che_gen['A+'].map(lambda x: pd.notnull(x))
    C_require = he_che_gen['A-'].map(lambda x: pd.notnull(x))
    Sa_require = he_che_gen['B'].map(lambda x: pd.notnull(x))
    Pr_require =he_che_gen['E+'].map(lambda x: pd.notnull(x))
    Pa_require =he_che_gen['E-'].map(lambda x: pd.notnull(x))
    J_require =he_che_gen['N'].map(lambda x: pd.notnull(x))
    Mp_require =he_che_gen['O'].map(lambda x: pd.notnull(x))
    J_require1=he_che_gen['K'].map(lambda x: pd.notnull(x))
    Mp_require1 =he_che_gen['Z'].map(lambda x: pd.notnull(x))
    some = he_che_gen[T_require |C_require| Sa_require | Pr_require | Pa_require|J_require|Mp_require|J_require1|Mp_require1]


    #删除DBID为空的对
    Db1_require = some['DB_ID1'].map(lambda x: x!='null')
    Db2_require = some['DB_ID2'].map(lambda x: x!='null')

    triples_che_gen = some[Db1_require & Db2_require]


    print("过滤合并完成")

    #提取出实体,并修改实体DB名称
    triples_che_gen ["DB_ID1"] = triples_che_gen ["DB_ID1"] .apply(lambda x :"C:"+str(x))
    triples_che_gen ["DB_ID2"] = triples_che_gen ["DB_ID2"] .apply(lambda x :"G:"+str(x))

    triples_che_gen.to_csv("triples/triples_che_gen_themes.tsv",sep="\t",header=True,index=False)
    print(triples_che_gen)

    che_gen_chemical= triples_che_gen[["Entity1",'DB_ID1']].drop_duplicates(['DB_ID1'])
    che_gen_gene = triples_che_gen[["Entity2",'DB_ID2']].drop_duplicates(['DB_ID2'])

    print(che_gen_chemical)

    che_gen_chemical.to_csv("original entitys/entity_cg_chemical.tsv",sep='\t',header=True,index=False)
    print("to_csv chemical")
    che_gen_gene.to_csv("original entitys/entity_cg_gene.tsv",sep='\t',header=True,index=False)
    print(che_gen_gene)
    print("to_csv chemical")

def relation_normalization():
#提取出三元组
    triple_che_gen= pd.read_csv("triples/triples_che_gen_themes.tsv",delimiter='\t')
    print(triple_che_gen.describe())

    print("最大值：",max(triple_che_gen[['N']].values))
    #45222
    #将关系进行归一化：
    #A+	A+.ind	A-	A-.ind	B	B.ind	E+	E+.ind	E-	E-.ind	E	E.ind	N	N.ind	O	O.ind	K	K.ind	Z	Z.ind
    #
    triple_che_gen ["A+"] = triple_che_gen ["A+"] .apply(lambda x : x/45222.0)
    triple_che_gen ["A-"] = triple_che_gen ["A-"] .apply(lambda x : x/45222.0)
    triple_che_gen ["B"] = triple_che_gen ["B"] .apply(lambda x : x/45222.0)
    triple_che_gen ["E+"] = triple_che_gen ["E+"] .apply(lambda x : x/45222.0)
    triple_che_gen ["E-"] = triple_che_gen ["E-"] .apply(lambda x : x/45222.0)
    triple_che_gen ["E"] = triple_che_gen ["E"] .apply(lambda x : x/45222.0)
    triple_che_gen ["N"] = triple_che_gen ["N"] .apply(lambda x : x/45222.0)
    triple_che_gen ["O"] = triple_che_gen ["O"] .apply(lambda x : x/45222.0)
    triple_che_gen ["K"] = triple_che_gen ["K"] .apply(lambda x : x/45222.0)
    triple_che_gen ["Z"] = triple_che_gen ["Z"] .apply(lambda x : x/45222.0)
    print("归一化完成")
    triples= []
    for index ,row  in triple_che_gen.iterrows():
        if index %50000 ==0:
            print("读入",index/50000,"行")
        entity_1 = row["DB_ID1"]
        entity_2 = row["DB_ID2"]
        A1 = row['A+']
        A2 = row['A-']
        B = row['B']
        E1 = row['E+']
        E2 = row['E-']
        E = row['E']
        N = row['N']
        O = row['O']
        K = row['K']
        Z = row['Z']
        if A1!=0.0:
            triples.append([entity_1,"A+",entity_2,A1])
        if A2 != 0.0:
            triples.append([entity_1, "A-", entity_2, A2])
        if B != 0.0:
            triples.append([entity_1, "B", entity_2, B])
        if E1 != 0.0:
            triples.append([entity_1, "E+", entity_2, E1])
        if E2 != 0.0:
            triples.append([entity_1, "E-", entity_2, E2])
        if E != 0.0:
            triples.append([entity_1, "E", entity_2, E])
        if N != 0.0:
            triples.append([entity_1, "N", entity_2, N])
        if O != 0.0:
            triples.append([entity_1, "O", entity_2, O])
        if K != 0.0:
            triples.append([entity_1, "K", entity_2, K])
        if Z != 0.0:
            triples.append([entity_1, "Z", entity_2, Z])
    print("read done")
    che_gen_triples = DataFrame(triples,columns=["che","rel","gen","score"])

    #去重
    che_gen_triples.drop_duplicates(["che","rel","gen","score"],inplace=True)
    final_che_gen_triples = che_gen_triples.sort_values('score', ascending=False).groupby(["ent1","rel","ent2"]).first().reset_index()
    final_che_gen_triples.to_csv("triples_che_rel_gen.tsv",sep='\t',header=False,index=False)
    print(che_gen_triples)

if __name__ == '__main__':
    part1_pat2_concat()
    relation_normalization()