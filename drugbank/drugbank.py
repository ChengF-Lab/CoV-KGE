import pandas as pd
pd.set_option("precision",10)
def dbgn():
    db = pd.read_csv("drugbank_all_triples.tsv",delimiter='\t',names=["db","rel","ent2"])
    db_ent = pd.read_csv("drugbank_entity.tsv",names=["db"],delimiter='\t')

    #把GNBR中有对应的drugbank id的所有meshid以及chebi id 实体找出来
    #mesh-drugbank的映射文件
    mesh_db = pd.read_csv("mesh_drugank.tsv",delimiter='\t',names = ["db","gnbrid"])

    # print(mesh_db["gnbrid"].value_counts())
    # print(mesh_db["db"].value_counts())
    GNBR_ent = pd.read_csv("../entity_GNBRid.tsv",delimiter='\t',names=["id","name","gnbrid"])
    mesh  = GNBR_ent[GNBR_ent["gnbrid"].str.contains("C:MESH:")]
    gnbr_mesh_db = mesh.merge(mesh_db,on="gnbrid",how='left')
    print(gnbr_mesh_db)
    in_gnbr_mesh_db = gnbr_mesh_db[~gnbr_mesh_db["db"].isnull()]
    print(in_gnbr_mesh_db)
    ss = in_gnbr_mesh_db[["gnbrid","db"]]
    ss.to_csv("gnbr_mesh_db.tsv",sep='\t',header=True,index=False)

    chebi  = GNBR_ent[GNBR_ent["gnbrid"].str.contains("C:CHEBI:")]
    #chebi-drugbank的映射文件
    chebi_db = pd.read_csv('../drug_finding_6_rel_filter_middle/drugbank_chebi.tsv',delimiter='\t')

    chebi_db["gnbrid"] = chebi_db["gnbrid"].apply(lambda x: "C:CHEBI:"+str(x))
    # chebi_db.rename({"chebi":"gnbrid","db":"drugbank"},inplace=True)
    print(chebi_db["gnbrid"].value_counts())
    print(chebi_db["db"].value_counts())
    print(chebi_db)
    gnbr_chebi_db= chebi.merge(chebi_db,on="gnbrid",how="left")
    in_gnbr_chebi_db = gnbr_chebi_db[~gnbr_chebi_db["db"].isnull()]
    qq = in_gnbr_chebi_db[["gnbrid","db"]]
    qq.to_csv("gnbr_chebi_db.tsv",sep='\t',header=True,index=False)

    print(qq["gnbrid"].value_counts())

    # part 2 把所有有对应 dbid 的实体合并起来，并去重
    mesh_in = pd.read_csv("gnbr_mesh_db.tsv",delimiter='\t')
    chebi_in = pd.read_csv("gnbr_chebi_db.tsv",delimiter='\t')
    db_in_gnbr= pd.concat([mesh_in,chebi_in],axis=0)
    print(db_in_gnbr)
    db_in_gnbr.drop_duplicates(["db"],keep="first",inplace=True)
    print(db_in_gnbr)
    db_in_gnbr.to_csv("gnbr_all_db.tsv",sep='\t',header=False,index=False)


    # part2 将DDI中的dbid 转换成gnbr的id：去掉不在gnbr中的ddi
    db_in_gnbr = pd.read_csv("gnbr_all_db.tsv",delimiter='\t',names=["gnbr","db"])
    db_in_gnbr["db"] = db_in_gnbr["db"].apply(lambda  x: "<http://bio2rdf.org/drugbank:"+x+">")
    print(db_in_gnbr)
    in_db = db_in_gnbr[["db"]].iloc[:,:].values
    db_triples= pd.read_csv("drugbank_all_triples.tsv",delimiter='\t',names=["ent1","rel","ent2"])

    ddi = db_triples[db_triples["rel"] == "ddi-interactor-in"]
    ddi = ddi[db_triples["ent1"].isin(in_db[:,0]) & ddi["ent2"].isin(in_db[:,0])]
    print(ddi)
    first = ddi.merge(db_in_gnbr,left_on="ent1",right_on="db",how = "left")
    print(first)
    first.to_csv("left.tsv",sep = '\t',header = True,index=False)
    left = first[["gnbr","rel","ent2"]]
    left.rename({"gnbr":"left"},inplace=True)

    right = left.merge(db_in_gnbr,left_on="ent2",right_on="db",how='left')
    print(right)
    right = right[["gnbr_x","rel","gnbr_y"]]
    right.to_csv("right.tsv",sep='\t',header=False,index=False)
    ss=right.drop_duplicates(subset=["gnbr_x","gnbr_y"])
    print(ss)

    #把属性三元组添加进来（GNBR中有drugbank id的药物的属性三元组找出来）
    other = db_triples[~db_triples["rel"] .isin(["ddi-interactor-in"]) ]
    othre_triple = other[other["ent1"].isin(in_db[:,0]) | other["ent2"].isin(in_db[:,0])]
    print(othre_triple)
    first = othre_triple.merge(db_in_gnbr,left_on="ent1",right_on="db",how = "left")
    aa = first[["gnbr",'rel',"ent2"]]
    aa.to_csv("shuxin_triples.tsv",sep = '\t',header=True,index=False)

    shuxin = pd.read_csv("shuxin_triples.tsv",delimiter='\t')
    print(shuxin)
    ent  = shuxin[["ent2"]].iloc[:,:].values
    entitys = list(set(ent[:,0]))
    entity = pd.DataFrame(entitys)
    ent.drop_duplicates(keep = "first",inplace=True)
    print(ent)
    #将属性三元组中的出现的关系提取出来
    print(entity)
    entity.to_csv("new_add_ents.tsv",sep='\t',header=False,index=False)
    relations = shuxin[["rel"]].iloc[:,:].values
    re =list(set( relations[:,0]))
    res = pd.DataFrame(re)
    print(res)
    res.to_csv("new_add_relations.tsv",sep = "\t",header = False,index = False)


    # 将所有三元组合并：ddi +属性三元组+gnbr三元组
    ddi = pd.read_csv("ddi_triples.tsv",delimiter='\t',names=["h","r","t"])
    print(ddi)
    shu = pd.read_csv("shuxin_triples.tsv",delimiter='\t',names=["h","r","t"])
    print(shu)
    gnbr = pd.read_csv("final_gnbr_triples.tsv",delimiter='\t',names=["h","t","r","s"])
    print(gnbr)
    gnbl  = gnbr[["h","r","t","s"]]
    print(gnbl)
    all_tri = pd.concat([ddi,shu,gnbl],axis=0)
    all_tri.to_csv("./dbgn/dbgn_triples.tsv",sep='\t',header=False,index=False)
    print(all_tri)

    # 所有实体合并：新增的只有属性中出现的实体
    new_add = pd.read_csv("new_add_ents.tsv",delimiter='\t',names=["h"])
    print(new_add)
    gn = pd.read_csv("final_gnbr_entities.tsv",delimiter='\t',names=["h"])
    print(gn)
    all_ent = pd.concat([new_add,gn],axis=0)
    all_ent.to_csv("./dbgn/dbgn_entities.tsv",sep='\t',header=False,index=True)
    print(all_ent)
    C = all_ent[all_ent["h"].str.contains("C:")]
    print(C)
    G= all_ent[all_ent["h"].str.contains("G:")]
    print(G)
    D= all_ent[all_ent["h"].str.contains("D:")]
    print(D)
    #改变id 顺序

    new_add = pd.read_csv("./dbgn/dbgn_entities.tsv",delimiter='\t',names=["id","h"])
    all_ent = new_add[["h","id"]]
    all_ent.to_csv("./dbgn/dbgn_entities.tsv",sep='\t',header=False,index=False)

    #所有关系合并
    relation = pd.read_csv("../gnbr/relations.tsv", delimiter='\t', names=["rel","id"])
    rel1 = relation[["rel"]]
    rel2 = pd.read_csv("new_add_relations.tsv", delimiter='\t', names=["rel"])
    all_rel = pd.concat([rel1, rel2], axis=0)
    all_rel.to_csv("./dbgn/dbgn_relations.tsv", sep='\t', header=False, index=True)

    #统计各类关系三元组数量
    # gnbr = pd.read_csv("final_0326_gnbr.tsv",delimiter='\t',names=["h","t","r","s"])
    # cg = gnbr[gnbr["h"].str.contains("C:") & gnbr["t"].str.contains("G:")]
    # print("cg",cg)
    # che_dis = gnbr[gnbr["h"].str.contains("C:") & gnbr["t"].str.contains("D:")]
    # print("che_dis",che_dis)
    # gen_gen = gnbr[gnbr["h"].str.contains("G:") & gnbr["t"].str.contains("G:")]
    # print("gen_gen",gen_gen)
    # gen_dis = gnbr[gnbr["h"].str.contains("G:") & gnbr["t"].str.contains("D:")]
    # print("gen_dis",gen_dis)


if __name__ == '__main__':
    dbgn()