import pandas as pd
def main():
    df=pd.read_excel("microsoft_paper_nasties.xlsx")
    #print(df.columns)
    df=df[["source","chain"]]
    split_source=df["source"].str.split(" ",expand=True)
    split_source.columns=["db","id"]
    split_source["db"]=split_source["db"].str[:-1]
    df["db"]=split_source["db"]
    df["id"]=split_source["id"]
    df=df.drop(columns="source")
    print(df)
    df.to_csv("capiti_E_splitdb.tsv",sep="\t")
    #df=df.rename({"Sequence Source ":"source","Chain":"chain"},axis=1)

if __name__ =="__main__":
    main()
