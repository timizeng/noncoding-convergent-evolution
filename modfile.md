.............通过4d位点获得的non.cons.mod, 均以Chr1为例..............
```perl
ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -1606820.522595
BACKGROUND: 0.276160 0.184605 0.149641 0.389594
RATE_MAT:
  -0.942714    0.216777    0.346913    0.379023
   0.324288   -1.285675    0.222231    0.739157
   0.640224    0.274157   -1.329930    0.415548
   0.268667    0.350242    0.159609   -0.778519
TREE: (((arabidopsis_thaliana:0.730447,theobroma_cacao:0.414062):0.0180751,(citrus_sinensis:0.276673,(swietenia_macrophylla:0.115554,(xylocarpus_granatum:0.0646164,(xylocarpus_rumphii:0.0782978,xylocarpus_moluccensis:0.0691567):0.0104676):0.0622751):0.13709):0.147878):0.0481315,(populus_trichocarpa:0.351873,(carallia_pectinifolia:0.170954,((bruguiera_gymnorhiza:0.0671421,bruguiera_sexangula:0.0812895):0.0486947,((ceriops_tagal:0.128577,(kandelia_obovata:0.0496971,kandelia_candel:0.0528404):0.0684067):0.0151461,(rhizophora_mangle:0.074894,(rhizophora_apiculata:0.0634833,(rhizophora_mucronata:0.063605,rhizophora_stylosa:0.0507747):0.0144585):0.0272553):0.0523516):0.0156791):0.0580341):0.219283):0.0481315);
```


.....通过estimate-rho,--rho 0.4 --target-coverage 0.125 --expected-length 20 计算的 non.cons.mod......
```perl
ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -17296270.035638
BACKGROUND: 0.276160 0.184605 0.149641 0.389594
RATE_MAT:
  -0.942714    0.216777    0.346913    0.379023
   0.324288   -1.285675    0.222231    0.739157
   0.640224    0.274157   -1.329930    0.415548
   0.268667    0.350242    0.159609   -0.778519
TREE: (((arabidopsis_thaliana:0.730447,theobroma_cacao:0.414062):0.0180751,(citrus_sinensis:0.276673,(swietenia_macrophylla:0.115554,(xylocarpus_granatum:0.0646164,(xylocarpus_rumphii:0.0782978,xylocarpus_moluccensis:0.0691567):0.0104676):0.0622751):0.13709):0.147878):0.0481315,(populus_trichocarpa:0.351873,(carallia_pectinifolia:0.170954,((bruguiera_gymnorhiza:0.0671421,bruguiera_sexangula:0.0812895):0.0486947,((ceriops_tagal:0.128577,(kandelia_obovata:0.0496971,kandelia_candel:0.0528404):0.0684067):0.0151461,(rhizophora_mangle:0.074894,(rhizophora_apiculata:0.0634833,(rhizophora_mucronata:0.063605,rhizophora_stylosa:0.0507747):0.0144585):0.0272553):0.0523516):0.0156791):0.0580341):0.219283):0.0481315);
```

.....通过estimate-rho,--rho 0.4 --target-coverage 0.125 --expected-length 20 计算的 cons.mod......
```perl
ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -17296270.035638
BACKGROUND: 0.276160 0.184605 0.149641 0.389594
RATE_MAT:
  -0.942714    0.216777    0.346913    0.379023
   0.324288   -1.285675    0.222231    0.739157
   0.640224    0.274157   -1.329930    0.415548
   0.268667    0.350242    0.159609   -0.778519
TREE: (((arabidopsis_thaliana:0.163551,theobroma_cacao:0.0927108):0.00404711,(citrus_sinensis:0.0619486,(swietenia_macrophylla:0.0258732,(xylocarpus_granatum:0.014468,(xylocarpus_rumphii:0.0175313,xylocarpus_moluccensis:0.0154846):0.00234375):0.0139437):0.0306952):0.0331107):0.0107769,(populus_trichocarpa:0.0787863,(carallia_pectinifolia:0.0382775,((bruguiera_gymnorhiza:0.0150335,bruguiera_sexangula:0.0182012):0.010903,((ceriops_tagal:0.0287891,(kandelia_obovata:0.0111275,kandelia_candel:0.0118313):0.0153166):0.0033913,(rhizophora_mangle:0.0167692,(rhizophora_apiculata:0.0142143,(rhizophora_mucronata:0.0142415,rhizophora_stylosa:0.0113687):0.00323734):0.00610261):0.0117218):0.00351064):0.0129942):0.0490987):0.0107769);
```


............通过滑窗计算的ave.noncons.mod, --gc 0.4 -target-coverage 0.125 --expected-length 20..........
```perl
ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
BACKGROUND: 0.300000 0.200000 0.200000 0.300000
RATE_MAT:
  -0.913147    0.172278    0.471411    0.269458
   0.258417   -1.129295    0.167011    0.703867
   0.707116    0.167011   -1.133529    0.259402
   0.269458    0.469245    0.172934   -0.911637
TREE: (((arabidopsis_thaliana:0.283559,theobroma_cacao:0.201499):0.033044,(citrus_sinensis:0.157516,(swietenia_macrophylla:0.090951,(xylocarpus_granatum:0.050449,(xylocarpus_rumphii:0.0557602,xylocarpus_moluccensis:0.0537152):0.0136772):0.0493806):0.075586):0.0687828):0.0228462,(populus_trichocarpa:0.186472,(carallia_pectinifolia:0.124836,((bruguiera_gymnorhiza:0.0587184,bruguiera_sexangula:0.0631935):0.0391498,((ceriops_tagal:0.0905301,(kandelia_obovata:0.0448993,kandelia_candel:0.0450836):0.0494373):0.0162342,(rhizophora_mangle:0.0721005,(rhizophora_apiculata:0.0494136,(rhizophora_mucronata:0.0498791,rhizophora_stylosa:0.0441812):0.0136202):0.0211542):0.0426692):0.0124291):0.0497753):0.0964312):0.0228462);
```

..........通过滑窗计算的ave.cons.mod, --gc 0.4 -target-coverage 0.125 --expected-length 20..............
```perl
ALPHABET: A C G T
ORDER: 0
SUBST_MOD: REV
BACKGROUND: 0.300000 0.200000 0.200000 0.300000
RATE_MAT:
  -0.913147    0.172278    0.471411    0.269458
   0.258417   -1.129295    0.167011    0.703867
   0.707116    0.167011   -1.133529    0.259402
   0.269458    0.469245    0.172934   -0.911637
TREE: (((arabidopsis_thaliana:0.0883812,theobroma_cacao:0.0623101):0.0102761,(citrus_sinensis:0.0487847,(swietenia_macrophylla:0.0279154,(xylocarpus_granatum:0.0154625,(xylocarpus_rumphii:0.0172362,xylocarpus_moluccensis:0.0165067):0.00412912):0.0152416):0.0233294):0.0214665):0.00713203,(populus_trichocarpa:0.0576809,(carallia_pectinifolia:0.0383625,((bruguiera_gymnorhiza:0.0180177,bruguiera_sexangula:0.0195072):0.0120448,((ceriops_tagal:0.0281078,(kandelia_obovata:0.0137497,kandelia_candel:0.0138919):0.015195):0.00495317,(rhizophora_mangle:0.022015,(rhizophora_apiculata:0.0152364,(rhizophora_mucronata:0.0153947,rhizophora_stylosa:0.0135783):0.00420926):0.00655642):0.0131142):0.00384291):0.0151912):0.0300964):0.00713203);

```

