import pandas as pd
pd.set_option("precision",10)
def triples():
    che_dis = pd.read_csv("triples_che_rel_dis", delimiter='\t', names=["ent1", 'rel', 'ent2', 'score'])

    che_gen = pd.read_csv("triples_che_rel_gen", delimiter='\t', names=["ent1", 'rel', 'ent2', 'score'])
    gen_dis = pd.read_csv("triples_gen_rel_dis.tsv", delimiter='\t', names=["ent1", 'rel', 'ent2', 'score'])
    gen_gen = pd.read_csv("triples_gen_rel_gen.tsv", delimiter='\t', names=["ent1", 'rel', 'ent2', 'score'])

    triples =pd.concat([che_dis,che_gen,gen_dis,gen_gen])
    print("done")
    triples.to_csv("triples.tsv", sep="\t", header=False, index=False)

    #得到所有实体
    # dis
    entity_cd_disease = pd.read_csv("original entitys/entity_cd_disease.tsv", delimiter="\t", names=None)
    print("done")
    print(entity_cd_disease)
    entity_gd_disease = pd.read_csv("original entitys/entity_gd_disease.tsv", delimiter="\t", names=None)
    print("done")
    print(entity_gd_disease)

    new_dis = pd.concat([entity_cd_disease, entity_gd_disease], axis=0)
    new_dis.rename(columns={'Entity2': 'entity', 'DB_ID2': 'GNBR_id'}, inplace=True)
    filter_disease = new_dis[["entity", 'GNBR_id']].drop_duplicates(['GNBR_id'])
    print(filter_disease)

    ######################
    # chemical
    entity_cd_chemical = pd.read_csv("original entitys/entity_cd_chemical.tsv", delimiter="\t", names=None)
    print("done")
    print(entity_cd_chemical)
    entity_cg_chemical = pd.read_csv("original entitys/entity_cg_chemical.tsv", delimiter="\t", names=None)
    print("done")
    print(entity_cg_chemical)

    new_che = pd.concat([entity_cd_chemical, entity_cg_chemical], axis=0)
    new_che.rename(columns={'Entity1': 'entity', 'DB_ID1': 'GNBR_id'}, inplace=True)
    filter_chemical = new_che[["entity", 'GNBR_id']].drop_duplicates(['GNBR_id'])
    print(filter_chemical)

    ##############################
    # gene
    entity_cg_gene = pd.read_csv("original entitys/entity_cg_gene.tsv", delimiter="\t")
    entity_cg_gene.rename(columns={'Entity2': "entity", 'DB_ID2': 'GNBR_id'}, inplace=True)
    print("done")
    print(entity_cg_gene)
    entity_gd_gene = pd.read_csv("original entitys/entity_gd_gene.tsv", delimiter="\t")
    entity_gd_gene.rename(columns={'Entity1': "entity", 'DB_ID1': 'GNBR_id'}, inplace=True)
    print("done")
    print(entity_gd_gene)

    entity_gg_gene1 = pd.read_csv("original entitys/entity_gg_gene1.tsv", delimiter="\t")
    entity_gg_gene1.rename(columns={'Entity1': "entity", 'DB_ID1': 'GNBR_id'}, inplace=True)
    print("done")
    print(entity_gg_gene1)
    entity_gg_gene2 = pd.read_csv("original entitys/entity_gg_gene2.tsv", delimiter="\t")
    entity_gg_gene2.rename(columns={'Entity2': "entity", 'DB_ID2': 'GNBR_id'}, inplace=True)
    print("done")
    print(entity_gg_gene2)
    new_gene = pd.concat([entity_cg_gene, entity_gd_gene, entity_gg_gene1, entity_gg_gene2], axis=0)
    filter_gene = new_gene[["entity", 'GNBR_id']].drop_duplicates(['GNBR_id'])
    print(filter_gene)
    entity_to_GNBRid = pd.concat([filter_disease, filter_chemical, filter_gene], axis=0)
    print("all_entities：", entity_to_GNBRid)
    entity_to_GNBRid.to_csv("entity_GNBRid.tsv", sep="\t", header=True, index=False)

    entities_id = entity_to_GNBRid["GNBR_id"]

    entities_id.to_csv("gnbr_entities.tsv", sep="\t", header=True, index=False)

def prepare():
    #去重特殊的chemical实体
    gnbr_tri = pd.read_csv("triples.tsv", delimiter='\t', names=["h", "r", "t", "score"])
    gnbr_ent1 = pd.read_csv("gnbr_entities.tsv", delimiter='\t', names=["gnbrid"])
    gnbr_ent1 = gnbr_ent1[["gnbrid","id"]]

    not_gnbrc = gnbr_ent1[~gnbr_ent1["gnbrid"].str.contains("C:")]
    gnbr_ent = gnbr_ent1[gnbr_ent1["gnbrid"].str.contains("C:")]

    c_ent = gnbr_ent[~(gnbr_ent["gnbrid"].str.contains("C:MESH") | gnbr_ent["gnbrid"].str.contains("C:CHEBI"))]
    no_c_ent = gnbr_ent[(gnbr_ent["gnbrid"].str.contains("C:MESH") | gnbr_ent["gnbrid"].str.contains("C:CHEBI"))]
    print(c_ent)
    cs = c_ent[["gnbrid"]].iloc[:, :].values

    c_ent["gnbrid"] = c_ent["gnbrid"].apply(lambda x: x.replace(":", ":MESH:"))
    print(c_ent)
    all_ent = pd.concat([no_c_ent, not_gnbrc, c_ent], axis=0)
    print(all_ent)
    all_ent.drop_duplicates(["gnbrid"], keep="first", inplace=True)
    all_ent = all_ent["gnbrid"]
    print(all_ent)
    all_ent.to_csv("final_gnbr_entities.tsv", sep='\t', header=False, index=False)
    # 把C开头的三元组提取出来
    # not_c = gnbr_tri[~(gnbr_tri["h"].str.contains("C:") | gnbr_tri["t"].str.contains("C:"))]

    not_ctri = gnbr_tri[~gnbr_tri["h"].isin(cs[:, 0])]
    print(not_ctri)
    c_tri = gnbr_tri[gnbr_tri["h"].isin(cs[:, 0])]
    print(c_tri)
    c_tri["h"] = c_tri["h"].apply(lambda x: x.replace(":", ":MESH:"))
    print(c_tri)
    all_tri = pd.concat([not_ctri, c_tri], axis=0)
    print(all_tri)
    final_tri = all_tri.sort_values(by=["score"], ascending=False).groupby(by=["h", "t", "r"]).first().reset_index()
    print(final_tri)
    final_tri.to_csv("final_gnbr_triples.tsv", sep='\t', header=False, index=False)


if __name__ == '__main__':
    triples()
    prepare()


