## Table of Contents

- [Users' Guide](#uguide)
  - [General usage](#general)
  - [Making predictions](#prediction)

## <a name="uguide"></a>Users' Guide

snpMAP is a single-nucletide-variant (SNP) interpretation program that works on synonymous-SNPs and non-synonymous-SNPs. It takes a VCF file in input and annotates it with the predictions.


### <a name="general"></a>Configuration

snpMAP requires a configuration file conformant with this standard:

```sh
{
    "classification": "binary",
    "use_polyphen": "True", # Set True to use polyphen features 
    "retrain": "False",
    "input":{
        "data_dir": "data",
        "dependencies_dir": "dependencies",
        "case_path": "/home/anej/repos/studies/CBM/P2/data/case_processed_v2.vcf",  # Input VCF file
        "clinvar_path": "/home/anej/repos/studies/CBM/P2/train_data/clinvar.vcf" # SNP clinvar annotations (for train)
    },
    "output":{
        "cache_dir": "cache",
        "results_dir": "results",
        "classification_dir": "classification"
    },
    "SeattleSeq":{ # Set the java classpaths for the 4 jars in the dependecies and the seattle138 dir. 
        "java_classpath": "./
        :/home/anej/repos/studies/CBM/P2/dependencies/seattle138/httpunit-1.7.jar
        :/home/anej/repos/studies/CBM/P2/dependencies/seattle138/js-1.6R5.jar
        :/home/anej/repos/studies/CBM/P2/dependencies/seattle138/nekohtml-0.9.5.jar
        :/home/anej/repos/studies/CBM/P2/dependencies/seattle138/xercesImpl-2.6.1.jar
        :/home/anej/repos/studies/CBM/P2/dependencies/seattle138/"
    }
}


```


### <a name="prediction"></a>Making predictions

Once the config.json file is configured, annotating your VCF file with predictions is as simple as:

```sh
python3 snpMAP.py     
```

If you want to train the model on your machine, please set the "retrain" option to True, otherwise pre-trained model will be used for prediction. This, however, is maily used for developmental purposes, thus, is not needed and will result the same model as provided in the repository.

The use_polyphen flag should be set to True by default and only disabled if the PolyPhen2 server is not responding (for more than 20-30 minutes). Disabling it allows for making the predictions without the PolyPhen2 features.

The default mode of the model is "binary". The use of "multiclass" is possible for the nsSNP's, but does not result in a substantially better predictor, and thus, the corresponding model is not included in the repository.