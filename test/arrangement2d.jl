using Plasm


# //////////////////////////////////////////////////////////////////////////////
function generate_debugging_data()

  V = [0.8006621680849682 0.7409967402555993 0.3948057913681639 0.1189697641436624 0.21328131471412914 0.25906885623179376 0.9734091678813379 0.6776474971408346 0.6552239722511268 0.7895551412151476 0.07319596480351402 0.35933396436113524 0.021531116592966734 0.2578385241362297 0.37560727873054234 0.8765667083626952 0.6598793633195865 0.41945419274361495 0.09766019790032936 0.690072369891066 0.8777510109488597 0.9048137742731802 1.1133520680739248 0.7113509999567786 0.6048061081055881 0.16545422275076432 0.8529063482691522 0.3433813528249137 0.18084934687309895 0.1149782249632958 0.6864348737800479 0.007131734104217944 0.5536230469073733 0.5348784783909657 0.2179521988107871 -0.010910260107321545 0.7858399613753464 -0.030134084026466446 0.29416188082519645 0.5365187159151856 0.4743787155887167 0.8572836991903803 0.858578379604646 -0.04035841291025663 0.43475395251838445 0.7713087372153011 0.040439958778906694 0.7942363621366603 0.5905187150956013 1.0279093669525765 0.4406693493289661 0.856155167437376 0.5082737093310316 1.005019495273478 0.1073100543084187 0.6747038503697353 0.8873344032156302 0.6575203665835831 0.21086811379333809 0.8163455087912086 0.3438558496003455 0.8230959103871229 0.6127945924642452 0.7700914978238589 0.4927325170436216 0.27204687117855075 0.7142332939535492 0.5861840112300303 0.7434689799556513 0.4156091384224042 0.6396629854121297 0.48088280768901714 0.5751956778288142 0.7982971202171488 0.931443249679495 0.2331489366195833 0.8896669153826959 0.8117082101746101 0.36451376777702205 0.844506619050661 0.09345364593888154 0.29041693315627787 0.6706620660623666 0.40065904529022434 0.5004061073475131 0.8814333745602825 0.2751923568808375 0.19857657230719256 0.13511089041912133 0.13267135118984902 0.3862438915468885 0.2677466438180891 0.7789227355827202 0.7329442205589303 0.2954294866486687 0.6858353682989332 0.3924316073619485 0.9712694170644279 0.31627677732431614 0.523700073900955 0.9176176882365185 0.0925119205762081 -0.04517938670986965 0.5857941730126826 0.06871205763991442 0.22479196712136026 0.3841647987735649 0.14416052507569038 0.9494811392348863 0.24612762839213198 0.09186665510529635 0.43276168257732656 0.1506492102621005 0.4098829813727677 0.5830647533967397 0.6076492357465908 0.21306611718927598 0.8239520303005459 0.941447127877251 0.33901014778908134 0.8665487115171054 0.8135006429193083 0.8505639218484867 0.1437090219112202 0.39944479338148076 0.1408527612177753 0.5911577211780278 0.40213988030198966 0.2961789345513923 0.15877167825768515 0.46180058807852614 0.4047400923247756 0.6833993012948375 0.1688443606844146 0.5102967655334837 0.6960001961423984 0.2441218884085389 0.4979141700893497 0.1729390451253628 0.6621625017472256 0.1815187317732517 0.10306440023429045 0.7738203756289382 0.47258022571920033 0.9501625441160261 0.8449799949341082 0.2872250444919869 0.462897071987759 0.5541236077216544 0.6751987067153113 0.4666093329319599 0.4030627294049207 0.6297455951678427 0.7691709837671804 0.3572325188191825 0.7485161182628005 0.12873417519208857 0.17866755600201842 0.0360199818795265 0.7847584114189955 0.7487451643484574 0.25857545272591853 0.5166001209420114 0.4168620849418322 0.9141922158052822 0.6809911915534989 0.23899019430087126 0.5241106996909549 0.5391408280994063 0.2733330631167543 0.20819330136596592 0.8911499939668492 0.9219870957379972 0.36369093520606377 0.344507950343265 0.2452125214500909 0.16759265899472153 0.8412900208302171 0.6261844229192014 0.17961712976946312 0.10966731642890745 0.8225775193958815 0.46339146311605683 0.7829842189906248 0.24100629067246238 0.23645391481053082 0.56437333978591 0.881522097979062 0.13049765281763884 0.8305030649236494 0.1093694390872606 1.0258711246084231 0.6409388563723578 0.15138684163290708 0.4595897739231332 0.21374671008656934 0.31294009972450454 -0.0086509249108353 0.13939833610136598 0.3929534062127958 1.1302155531237368 0.4624315407876849 0.1795130277268701 0.9520945173592084 0.2048277817717032 0.2565993499701964 0.32852710337725816 0.38812789539289916 0.14035004939189927 0.8381267019206902 1.0115226402452253 0.6957673216767901 0.8125236735403396 0.1439741212714219 0.3657245830853695 0.8728721778881663 0.6148053681634134 1.1037233254619203 0.028603254670291733 0.6426164477653534 0.49230768484380016 0.2053443034238757 0.4772522951583003 0.12515621972148006 0.40805831394098235 0.8188433165318122 -0.014287644320647247 0.19036145635042545 0.082494334164912 0.6931182916542082 0.2588265782111159 0.5405457692235809 0.10195818915104085 0.023020071888582638 0.1822237214028031 0.9067616163874627 0.3298893426853525 0.7865872044021887 0.5728910516917506 -0.028125570702980544 0.3842547510637877 0.4611800402564069 0.3615659980119782 0.8288835084982249 0.7450680114671192 0.6703311454768539 0.5676318853065477 0.7458839566665351 0.3055739799922796 0.8937105489183593 0.056234486730007055 0.7557626724210095 0.32116657865803766 0.6704077936186972 0.49107497166486547 0.8287901330009412 0.7350406739238327 0.6209209907083433 0.8523931637332217 0.3249845472438296 0.9112145667582812 0.7587599780144743 0.5082273883913955 0.1683977530610847 0.28652908047470166 0.2003555350067954 0.8950242645855054 0.15185541586404722 0.12447064590432605 0.3635249659881429 0.8877169542604544 0.9573942957293297 0.6711633623914425 0.24694583731425285 0.03425940893387316 0.2459790821607312 0.8297578332272924 0.6048068091085136 0.23036070391403707 0.11580977404980791 -0.07105303900570746 0.5030725169594104 0.5872722697546383 0.6572364285827075 0.8091912248258345 0.8993808304184027 0.2924776017225664 0.7362916956905408 -0.07907436972632059 0.47746671317212075 0.22773820785476406 0.7925870981129921 0.6350467475778877 0.5644894126561988 0.7536030067869364 0.9194469376783901 0.6062128570416837 0.3572765638113934 0.31212588508764183 0.4660738443678369 1.0261496388500153 -0.05775099332887876 0.9033460760715273 0.48017763405713343 0.1652748327850598 0.5116149711324721 0.8367199242297878 0.2150391782369812 0.8692056867643613 0.8694864977676835 0.7749859904922172 0.6227423413947888 0.24394675372287228 0.8374135550379259 0.03447187623637667 0.7440284020270054 0.12590847768734237 0.8537630831757576 0.8874996359729143 0.19015198604534733 0.1439244798450108 1.036082019214568 0.7152551971219667 0.1391498346768403 0.6772545849603143 0.7937525356103263 0.935668120255525 0.3117578770912002 0.6458114477640124 0.4118098091840013 1.084689180888004 0.5850176192172838 0.7133807483260591 0.3653761411335219 0.6672437285400911 0.2228202389735861 0.5038760651866558 0.8885723238402821 0.1674387763292498 0.003017107913769256 0.09921283899656112 0.2909458674114092 0.0318888681990233 0.6568671579879874 0.9552786332717002 0.7087741278716014 0.9308162899366048 0.8539293308448693 0.49063082310065964 0.6101236277352197 0.06384226648353916 0.6060110650449927 0.3840564023081843 0.674534302276053 -0.11459855506650676 0.7463915614389804 0.8206901460139762 0.989997616598154 0.04596360748479087 0.5526657302843397 0.6527014822469823 0.25658060276461325 0.27072774933422616 0.11384817049761609 0.0960628626809632 0.8869171955832497 0.014923945608130085 0.0805211345579252 0.27069907863457054 0.9455448103916976 0.7765412421803493 0.2565569273821492 0.7090708637192334 0.2713555672534712 0.15188207971997408 0.578497177216539 0.7772523471077984 0.27952049177963206 0.2997478620670896 0.591505239569743 0.5328222592265207 0.6444502764478107 0.6355916286695921 0.3502256014893264 0.6836834541835645 0.25383483128837236 0.6312340102144656 0.67659975319196 0.9403407514822218 0.10863144378222822 0.2556026848472736 0.06407126876011152 0.3231267845370121 0.8790047550280033 0.5959299676275506 0.6004965226569688 0.15489005348298088 0.37871464313239567 0.8206071242124976 0.57457351648082 0.7703526774537088 0.7973208818549706 0.5461867109540877 0.4548580549708242 0.4058684119095108 0.06716876791219303 0.9701752618593955 0.9327268665127122 0.7553298449006789 0.9457501954341441 0.2848335449953781 0.3371212031276219 0.07433617673717159 0.2700263059575981 0.18118187208557168 1.0663977561276479 0.8664532158410436 0.2854034897593706 0.2707232644376125 0.49738059033711707 0.7398029976347613 0.8713314855700577 0.8469459729990432 0.6926259949659953 0.7852236135065531 -0.032505886804789635 0.8985120338171041 0.4174335194287532 0.4648970400001607 0.16487936522978192 0.59936128497454 0.04829664007945374 0.18881522979272275 0.8622796994528801 -0.006646417068877364 0.18550856195084198 0.9115068667416455 0.11197412749101032 0.4340900467957196 0.5089299734642215 0.42552040747404973 0.8354612550960143 0.8660205862039811 0.23322852549204418 0.6421642106758304 0.8018196592250084 -0.023264130455162894 0.6829906277825554 0.3287247691155617 0.8448519451268295 0.44183569071849754 0.7362995581453884 0.47977693698243684 0.9481803084510914 0.12918907159794 0.5875853489086772 0.8873610101684839 0.682641273580836 0.03703671904871225 0.8238981022387246 0.32255036702440293 0.9947820790916264 0.6138129234331179 0.8516689923654511 0.5755833541249968 0.26194062646389626 0.4786321543976675 0.46349801492840514 0.9527245905463476 0.49044503380186677 0.4106780557342861 0.7175416031853491 0.5968279658743908 0.654659010178085 0.9023502850106812 0.4021796064889345 0.6405999196501933 0.4731150812228827 0.19722517828414832 0.816313163366912 0.2895296321782489 0.41964117321282485 0.569988932657336 0.709246175177472 0.29944556509947107 0.8861592891465271 0.16178904441403347 0.42348555097540524 0.3470342691499709 0.2326824066128496 0.13746309671325613 0.34587348296135473 0.7182060433086408 0.7999601314107172 0.43565813237467654 0.5538332780316355 0.355278828688802 1.017298071059703 0.6410054912791959 0.8153411799059547 0.853546431018357 -0.013120482797944644 0.07203810820777606 0.581238946318887 0.10816732218525113 0.144210887205009 0.3628249778726068 -0.009967381897236172 0.7204670965571982 0.08271608213035377 -0.05372998523283823 0.3857181512045396 0.19463172454237535 0.4044727244106917 0.6502280086280807 0.5162455644185817 0.1997330864410969 1.0137206743566762 0.6709435790073329 0.21753791476372059 0.9691652771146967 0.7832896375884144 0.7516843303164168 -0.06603886905564499 0.58385852048839 0.35083053387680313 0.30531928653405094 0.42926572077455283 0.16635325438727921 0.41995484698247654 0.5028107175302315 0.5637202414119619 0.5487827718073875 0.06507000131481015 0.5190283349626784 0.7509559313519504 0.20137644579948163 0.45790456926781387 0.010711653221795181 0.8753149226378953 0.22918668650397814 0.14283162138599412 1.0639882564934784 0.3698306787905016 0.9273184858451401 0.9255105291233535 0.4634420912329171 0.5379049195203194 0.4032705517307219 0.6516190589823199 0.4551035290502738 0.6399002430750607 0.5449353778485573 0.3944035101963267 0.5147406749239447 0.7102695421684632 0.12033873752679632 0.49544599302862913 0.17539332577821437 0.6800644591561559 0.7010357230519667 0.16552863609911636 0.48361575265798107 0.7225858509783796 1.0149266937644672 0.7214635184183275 0.21683007467865306 0.3756572846212979 0.21403226237838105 0.09085253421267356 0.34626218631955086 0.8071185807852241 0.7186105804539917 0.33807323541970263 0.43484158011358576 0.09290544655677155 0.20727173555010803 0.801727297799259 0.5659811436147234 0.08430660337992688 0.022745449512923628 0.5734465437553478 0.6310480373416077 0.828273282096098 0.25475629599366 0.29867693405137674 0.7994722152195571 0.88888051729988 0.05001731761677672 0.8091400794451534 0.13858513646739 0.9459233225175764 0.5379245055564862 0.030160097311569295 0.5017232319051744 0.1662332401181015 0.32397303929891097 0.16080544905965155 0.1552078972419803 0.7043538204403983 0.8614222651917043 0.20078539802606654 0.20396815235728855 0.7750726103624925 0.42316312369278714 0.48532414262543383 0.2726767189825054 0.0816817072530055 0.0020920693369347737 0.9447540421464148 0.8576485046814873 0.688277026308425 0.7931196426931821 0.34211423859880136 0.23688028528821037 1.108431375677567 0.6587179997381618 0.7609010845595623 0.2569819984268577 0.3471180744723476 0.5535564152504364 0.21042851294227746 0.4598944028928175 0.3184476347967054 0.37977056571957996 0.6318908829706084 0.07465225577951926 0.07929167499254475 -0.038848056990723 0.6259277685725886 0.2591967452055672 0.7701444878909149 0.2229049354181631 0.01702706418320359 0.21605932442616213 0.8800280121307923 0.5741663996564025 0.8302002266064441 0.7831528092218581 0.06411516563576586 0.174790764523999 0.34153144189781653 0.16789005478326166 0.6017802258351108 0.7658876634319848 0.7031752369830349 0.7169344934963644 0.9434022139015178 0.4828481000251213 0.9080081028885284 0.3762898640825877 0.7906331397302896 0.3081924571318478 0.5558689319495344 0.47809432555654685 1.0025913769785182 0.4291323892923486 0.8296157415956187 0.9970501732581999 0.3561834322716674 0.9579629005025057 0.9482726536618653 0.8220697822684635 0.09862799791881274 0.26526942991773245 0.04045138870465423 0.8424671765315919 0.19818865981988215 0.024932237485693128 0.5131798141302256 0.6478994458247984 0.9334759849165214 0.666097870452056 0.22665197330868864 0.011040570234712143 0.23517571603825757 0.7536854495881466 0.75676787630053 0.21104710349248224 0.4086701882563675 0.1116542040561524 0.4596887512832979 0.9442309580654725 0.6320312498663769 0.8829490522972102 0.9687984061383904 0.10364758947655862 0.5865711782835747 0.2267472704039682 0.6073706051094897 0.4698442548380479 0.767144671603347 0.6792883125264667 0.7166823174501269 0.668297150643544 0.9392261905709122 0.609554973654599 0.328400874872805 0.4686681201482912 0.4989929626875924 0.7364288060160536 0.2226545453263918 0.7649519000467582 0.5205782175971004 0.11168248768990777 0.7739202432184581 1.0307126813975984 -0.017001486828491633 0.8300715190575961 0.7423142443938946 0.7919079330249801 0.538871003612106 0.06797336580956918 1.024837263712547 0.07355282608494966 0.5492526097691686 -0.09358695669290143 1.1150665226520828 0.8621683897467808 0.286943810706648 0.1073277148543265 0.7631946817735633 0.842939277990103 0.16068707738595053 0.6420435071173541 0.8451326690623169 0.848007658460044 0.2781820420730621 0.628790556868484 0.26732610477691393 0.8807626521343883 0.493042149549332 0.6066462058625092 0.38705555450133855 0.651035000636062 0.08570901074928143 0.3172680437802102 1.0436978682724458 -0.009921158737488833 0.06022930271068505 0.07130697026329076 0.19476741940205944 -0.006811251311319034 0.5595614214110142 0.9122431353064769 0.6345162164960877 0.9154795567284 0.871195913078063 0.7775918887842256 0.5437752060928235 0.09354540158955167 0.5122368896649724 0.6704469962724633 0.5015164418376608 0.14897271722035663 0.6323153400942143 0.7766519570816084 0.9847241278004012 -0.01891324843148621 0.5284876726773441 0.4925892866294105 0.33545861743079664 -0.0665638521987233 0.18079210121999856 0.0720797535633324 0.8164013589110451 -0.011034205699115552 -0.056156345578955436 0.45433996297881973 0.6435180406934901 0.47290843272191907 0.22445465897279526 0.7850125507679934 0.19982566927446815 -0.12529653061967322 0.6218927270988357 0.6513730530524683 0.3769020357523008 0.2347929378704484 0.3390114903181821 0.4813559863494265 0.7251199039333053 0.6658133547380506 0.1731779630959857 0.7434251700017306 0.10025086905587956 0.3555336248930612 0.5504429934401774 0.9081033357037487 0.11296146354521852 0.4194841564668999 0.12151929911985035 0.37503623425408844 1.0309441709506726 0.3852856364569969 0.5896059533904451 0.14559992134873026 0.16913209211598806 0.45845941572937765 0.9301637405924595; 0.9906836187663494 0.4120141109533818 0.1975406076447207 0.4999968889577415 0.846614354166266 0.6746625606555604 0.4444799647838364 0.29085734784705636 0.35490982591025766 0.2030697160406589 -0.09061557386487405 0.1414302206142438 0.21184100578864884 0.5747459499081751 0.4913985252893004 0.25039485716171966 0.13713589391960096 0.28055096026555865 0.9777894462837722 0.7090777520624281 0.9633288551692659 0.4218585319694082 0.12357781261311729 0.6090165016028964 0.09904771935338477 0.8238269090454593 0.7966904920365253 0.4104929750310523 0.6044322951399935 1.058653587763945 0.47539262005142885 0.900411549288661 0.08119611082360173 0.41272217432721425 0.6182352370177226 0.725645799833362 0.03382805731887068 0.9136465482472301 0.17815584457785102 0.09210337676290722 0.939743329651681 0.5351051341920547 0.8190138503486531 0.48773304109016663 0.27912554304989645 0.5982054848434817 0.7162536397562864 0.7485777379628966 1.0179743546585758 0.18235331907203364 0.42952867499547487 0.6412858920254122 0.10785324249947409 0.695144126942892 0.2688390225106439 0.09757324237192744 0.6265134357698201 0.0363694452381532 0.5520871883686383 0.9586061713531494 0.5193148914229624 -0.07191212278694403 0.08453984206674381 0.6334215652941179 0.7970996161213031 0.7398675217831351 0.7002863434519752 -0.027295744728543264 0.25491066995247036 0.9645820187147434 0.11137928521045751 0.4434523798797525 0.19575843459375333 0.6393602432211689 0.516905551720799 0.3594548842701849 0.8127525894350529 0.43984671684431875 0.4492351666563527 0.035702108963044085 0.43103849784108555 0.08246383052496131 0.5896540196641499 0.44331596952852054 0.35983242125422193 0.4541956707019819 0.4656326884532198 0.12424280525121491 0.15639885575023688 0.7586774091009897 1.0764004692959963 0.5406392636137758 -0.048897676045875235 0.5057301847306503 0.3598977678895471 0.010653521780280645 0.4019846232555271 0.3775323180659986 0.8313947322429023 0.9360246461347542 0.6557243339936414 0.8059690429547466 -0.011092495141713787 0.7953694122991986 0.08247751539494347 0.5900812187099592 0.9210828946728018 0.2129890699565996 0.48095810540545725 0.17081506917630457 0.9268494797837451 0.6559326926615723 0.8673550453318439 0.5492427408273256 0.4823643332428073 0.7538055694309346 0.8192813433583477 0.3440410901148345 0.09056659778638187 0.059057778756883915 0.314174691148849 0.7558092798169715 0.2668976969926198 0.51802724171018 1.0390450475046786 0.4414396411477047 0.41802063827604596 0.8892609914213275 0.7465316009138283 0.20187953300350767 0.4661902284796524 0.666859521583431 0.11495516112807247 0.9442421956372581 0.752348122406737 0.3145817255854754 0.7376476763552121 0.553424836674968 -0.14314671875317692 1.0410227179914118 0.9345178468199974 0.4268542711976102 0.4961413893011673 0.3703283092170201 0.6561216100031346 0.7315968963370918 0.14515403790163697 -0.07122387405218887 0.8855286313810371 0.4186480715522905 1.0450544683568825 0.21849868330231909 0.2610377938423699 0.21835813035030568 0.6114037281994266 0.5983927327007749 0.5913424254821927 0.6493402062764441 0.7004910407604341 0.4892046754564968 0.6388364273529512 0.12034185924163758 1.052934179974555 0.2109404872751073 0.9026355320707642 0.20041206953323615 0.7554799113143287 0.6746837814077187 0.2609883659408036 -0.04037451353078994 0.40381466550037043 0.3473770826282099 0.3193202606177419 0.6813084750481508 0.6459835766011994 0.9051318117845901 0.1738232741399168 -0.11280230897449905 0.3266768953290148 0.45892786538166586 0.2663103034748445 0.5898670100343448 0.11443390476025778 0.3295185153444249 0.7607594196835671 -0.06767862382169154 0.9340610138820928 0.38815050942882195 0.9030207629313873 0.12005133825633674 0.29989876061972537 0.09902561193123822 0.486056937391127 0.10640400766267474 0.8520896462473551 0.12664817759455602 0.5088103482157462 0.8597232355877412 0.22559319314823045 -0.04168577737693219 0.5832704454495992 0.2627073535152433 1.0447664976354454 0.6752011055191268 0.5518425826540809 0.7189359547017461 0.26765260630648924 0.6264688018474012 0.41628996643328847 0.2718112352074397 0.33117126572074257 0.2225489015179584 0.9791863433124737 0.9334123533454195 0.7089088800164951 -0.029518370453561343 0.10466464668514378 0.3501185462472858 0.14542898128318432 0.5994679239032286 0.48364311952881234 0.3919407509932555 0.7542742910118477 0.2114603198198592 0.13365314021876507 0.7293342900454803 0.5361275670019782 0.8663634815995372 0.3618923053659934 0.5423040803099048 0.5968608753939296 0.8912433840012857 0.5442629552387083 0.4800232757816286 0.45984089539086537 0.8335283211085812 0.913826780212621 0.42639913658726136 0.7831186451065503 0.6644744617821448 0.26361531864196125 0.9996012742371188 0.5911067803891918 0.24035005682831512 0.18603365562881982 0.6801504827551285 0.38026759020749373 0.945238694972456 0.7588932512835855 0.4627356215047949 0.9825056426429796 0.20602403070569864 0.36217950788108555 0.1702855285074919 0.4668884347157091 0.7408066018444092 0.6000455918587956 0.36101369039383624 0.1963288572061367 0.8710489113006826 0.9476418987141808 0.06728901538031062 0.6090193945593761 0.07480649448354806 0.09824338164997698 0.4356439418302218 0.2129218107375723 0.7143342416654526 0.4436118026052254 0.5360580309351678 0.9081422704570088 0.8764618691702408 0.9323133629853007 0.9322769421217989 0.7939442414181529 0.1463662503740547 0.44737902927456075 0.5602367790763433 0.3637246828321906 0.13369162207880597 0.6055985402895075 0.7666413717451445 0.424648846357549 1.1260452560137644 0.42172223235944206 0.8781085643939074 0.5546521161321397 0.22660392113281314 0.33882335239431727 0.41264730009158074 0.2612342525879349 0.5014618063109912 0.7712078747164529 0.7756532326897349 0.7676482782183793 0.7583398511392216 0.7253370123785523 0.994595823619469 0.3105107168365305 0.9490598894282161 0.7196859202562566 0.24141901437933083 0.9090542943450097 0.7604459002427473 0.2532905292845087 0.7383811320982784 0.6225556954321425 0.11102467743997854 0.5605991875961799 0.5699300374365109 0.08918971458124693 -0.006613518459344103 -0.06301600726151078 0.0993804760419375 0.6298303963259441 0.590688818267922 0.2603428500328235 0.15139429817049493 0.32397476849399265 0.8194775580871136 0.7001551193102031 0.864357026822445 1.0872165968886072 0.5446050578232284 0.8157227418344986 0.37056653195849026 0.7364166963437724 0.692968843249619 0.02939352044103717 0.32331738701394386 0.1635786852191261 0.9533914013182246 0.8145783373886596 0.4764601585932697 0.29563802564455866 0.5641442040445687 0.8721184152792499 0.567022853922618 0.2346296212123895 0.26897095744048016 0.6647317954077213 0.20785077277772204 0.40107035010202047 0.7607959291840515 0.8182578395750014 0.6244865333175863 0.288673913067069 0.7650431426083829 0.46303190959318075 0.5028219391711359 0.4606738643735905 0.6951139360819577 0.4446569233388506 0.906609607676092 0.35530677376459724 0.5148795671989693 0.6908848910611816 0.8031528650472712 0.14457486337311592 0.6112564348047972 0.0315769827212355 0.2620862839208826 1.0048597372027643 0.8187516001111762 0.2791059013763444 0.38882264959388946 0.21068383753682957 0.3757973188943235 0.6277355344125747 0.28775044727863663 0.8137907319834616 0.6980289888385309 1.0181335100174833 0.6966992792846648 0.206449500741066 0.6255485570110914 0.6966919033902262 0.5039539008660593 0.5623966830560004 0.559527402385258 0.7678017932237398 0.6976519414463678 0.6077556807341106 0.6873176822164925 0.9527862831692002 -0.018023591714267196 0.23164804582984358 0.8988245908406041 0.44691607482444046 0.886778699817472 0.9536153808182765 0.08091901087036638 0.3956899614579644 0.8308075408339983 0.48975746329718883 0.8252684387232502 0.0857836234825589 0.5385365267359226 0.18814744880020112 0.5224942449769515 0.9524338425177813 0.5195914739719555 0.02240789404209055 0.41859370676567065 0.7397696379076991 0.7106975028675303 0.54572366690297 0.48378562030086486 0.5428808953253244 -0.17414556482002427 0.0956721589631263 0.2375971411732151 0.4274873673290345 0.907277520940976 0.5854674692938379 0.1691217428019413 -0.0428466384877164 0.37847431515853064 0.6934231321574955 0.6683199720056325 0.7804981669670171 0.3162640035335588 0.08267500141818096 0.6711289199544003 -0.040391747751351376 0.7737195775402539 0.7210583329476526 0.48525286523300715 0.6903270328826691 0.8326102944089755 0.4610393051334095 1.0281175347951907 0.16308081839798766 0.4778552114704473 0.5123945121163666 0.6201911550853147 0.13587243822430467 0.9493922221419806 0.44881324599939965 0.37821464963129203 0.7812792860498793 0.4582058876602043 0.8057964015092588 0.30536540237317483 0.05422863779770905 0.3471935121038843 0.5647971118467188 0.8522233109084292 0.804923090671005 0.48518782868291344 0.581024692729038 0.8569013790276764 0.271062558948788 0.47986854371719395 0.3318652186023405 0.42837071035024554 0.7043649593734438 0.018052748474935897 0.5448918292064885 0.9723698916135579 0.4451835842709458 0.15556632853119495 -0.037697436935802346 0.7750629714657641 0.7761977288384028 0.712749592614553 0.8973596082370708 0.18132041060414208 0.3553111481964253 0.9937524006436551 0.0713441010727219 0.4404768808406916 -0.15257724072068982 0.6557539624001811 0.6496520651159843 0.48370723553525513 0.7042801963219485 0.4515653390785366 0.4011503473070169 0.15090458353197128 0.2543596645626826 0.287264518792618 0.5865388518981735 0.4859683045649152 0.15957368612959905 0.6759154932082412 0.42221697071239217 0.3203590019680361 0.07868272059377839 0.7749469878226006 0.7253572028394162 0.4242880744031392 0.11073259414849977 0.35008252303168597 0.5305997482859531 0.2644574340394247 0.590863827123002 0.49956118482881395 1.073594868123479 0.9183570881509707 0.5826413311483437 1.1542503324349902 0.14869290736084811 0.7482386250076534 0.40624164515847877 0.6010480152965009 0.8343484107192269 0.36500250465751255 0.36148521595931005 0.08766970307476424 0.9000991101017063 0.45249053551928786 0.9376872391227224 0.20063725300620608 0.7369125928735965 0.8088004164813362 0.6797120150068781 0.483904681300686 -0.06192615611638641 0.12692631319091072 0.27887390765898523 0.5186297309636236 0.3181108401923256 0.6306229453297937 0.6712673877813699 0.49823278948226374 0.3395383970513949 0.9080516819034102 1.0293329319225917 0.4914275552016132 0.23399712950780693 0.758320656998449 0.28101565866895556 1.0246364391111475 0.6124182138653107 0.30567636691383127 0.5550865904543866 0.5683819080089015 0.22409214701122643 0.9288703883857943 0.9815817335069056 0.3513650238842248 0.183131507713939 0.24678080535919725 0.6956619450693575 0.7395699181516104 -0.052854336852351996 0.15637469796300296 0.7956916332539259 0.3359717969278539 0.7625598050388734 0.1377811545308731 0.11508102588917846 0.3179724013578707 0.5380360726835342 0.5213848427328659 0.6188492591676792 0.7915082683407523 0.6586125568729597 0.7902180663920609 0.5277861550490888 0.2842764861112259 0.9290039194361221 -0.05290225727061479 1.0116712710553761 0.049672385086449944 0.952460250746886 0.589588855284117 0.14642245998234876 0.2876811045474946 0.2745017697462971 0.6044672509751283 0.42923767285232883 0.6538037313435912 0.8620418329942949 1.0887722827504325 0.39710614561246493 0.12809918785844732 0.4598275509058106 0.3059425526410845 0.22066568874804346 0.4407207061678091 0.33516327512252797 0.1592175142000074 0.6724047922185884 0.1206223962556487 1.0004160284402157 0.14800144516444297 0.9161114579775281 0.24677756624810415 0.2354623245125314 0.10117279917150084 0.449646098928296 0.0043955547544397774 0.7063076334129885 0.16031146184954584 0.607234348517808 0.6302921852995728 -0.011785118616270282 0.2211330442090474 0.6076750269336272 0.17716814826216026 0.9437874204544792 0.7103524848110522 0.5185713279109404 0.975241850679515 0.46987058219791145 0.6983054840619087 0.5604054653875568 0.08983387042606598 0.5407133158127896 0.1475956772566048 0.807121630407045 0.7998863413438075 0.8305718076533262 0.34327500799050736 -0.07309819027626621 0.504789302246643 0.28326138720925736 0.7855621260141088 0.33295336301424483 0.5863980053054512 0.7413375947889986 0.051740967651449796 0.29398770834438176 0.4398875352667754 0.5562945505243037 1.0179128529420303 0.3701994223303518 0.5086228378589672 0.7470675229273444 0.7089577536981579 0.521900166392126 0.4977025321452363 0.5982843680891731 0.9057341686019492 0.9829902820084706 0.2665569733056824 0.9248760339719327 0.3665347383197324 0.4919793316548066 0.9464842779743281 0.5629416164738859 0.4456180275328268 0.1809244890864545 0.7998923885923369 0.5172269551868242 0.8163297013577926 0.9897401220842834 0.36985239198455216 0.9361023905353106 0.444352288536344 0.1991663846725113 0.3194599237268882 0.7129898350108248 0.8218199125317728 0.7226450936467844 0.43827470669445756 0.13775474991171632 0.863310101404657 0.8979330426401793 0.2832730022274895 0.43879925464933817 0.08856371608952854 0.2548799557442059 0.5254263821064864 -0.039614231289644336 0.9504776220771913 0.48122676073127474 0.4646381849658157 0.6562482373449992 0.7259522659908555 0.7011263476744471 1.0303557179788145 0.6683687305175742 -0.04064865626129191 0.2997134201314373 0.6074138097751439 0.3712227079515241 -0.03119299926893522 0.719917693439272 0.7989467958427584 0.5917088052416188 0.8696550247807993 0.4357067945812443 0.749162683963624 0.8172239282736663 0.08291562710115097 0.3210352893865433 0.4443585347186813 0.07818619148530286 0.6350573977903442 0.3785858967454109 0.7072612976204128 0.9209893195671013 0.8612877772229521 0.6737693503696496 0.8490841464086962 -0.006483261415346164 0.8068265530356564 0.9521224175253155 0.2654903940089966 0.7049022697608021 0.47823897859682885 0.44413816875560386 0.7579601657699635 0.49038130335238467 0.23083442815628888 0.7708176189775989 0.6789245382575498 0.22284549385982577 0.31528794358348455 0.0811396311863871 -0.015539181854570976 0.8054892752891781 0.26157510431521447 0.3783628680425049 0.09973092797958893 0.6069592633939183 0.9194822877517039 0.9381247731448019 1.1023165995922812 0.8370184420503614 0.6042032279846458 0.9823753964951842 0.11279635776282709 0.8267492344506067 0.7956467810911325 0.18295367731452772 0.4497038693721549 -0.09347199200282927 0.9309883605778798 0.8032655909920959 0.4941981036217627 0.21055585005759112 0.8199594012487886 0.5990251106888677 0.6437911102167287 0.3766245248833424 0.33486016605808405 0.5598140885362325 0.19531558686148784 0.5322674560498827 0.771999626712778 1.058248902165611 0.84269370907373 0.485906732480323 0.7643449898858817 0.47750822253266123 0.26318830009507144 0.7877215126885929 0.7289646281907334 0.4674284067727846 0.8033165404776426 0.1642639333046269 0.7436815760601235 0.6402560947557909 0.5990011032109024 0.058188137993115015 0.6062603294279235 0.1258958861487101 0.16760077161253875 0.9912577009899081 1.115126623233091 0.13441661188845622 0.499652532822381 0.16174628606761404 0.5428723347312981 0.6777328342204426 0.26192848395805834 0.7827189848631182 0.6465812075315032 0.8517634653146513 0.5330720782841497 -0.008532756964207036 0.5504159545904924 0.36836940595666723 0.6074684319903033 0.625584407256846 0.5621411985288581 1.0451479709978362 0.7521294001826471 0.5986735929780979 0.8895360499504718 0.9239975861022138 0.03176350418024285 0.4998399477101858 1.0947717866521929 0.26415291676898384 1.0306661417931058 1.0231788406677447 0.4342942017545044 0.1764544732712688 0.7677037870618997 0.5774100856746247 1.0031230230878485 0.23419413876540002 0.8107943752523142 0.1166060103460208 0.6440993972090835];
  EV = [[1, 401], [2, 402], [3, 403], [4, 404], [5, 405], [6, 406], [7, 407], [8, 408], [9, 409], [10, 410], [11, 411], [12, 412], [13, 413], [14, 414], [15, 415], [16, 416], [17, 417], [18, 418], [19, 419], [20, 420], [21, 421], [22, 422], [23, 423], [24, 424], [25, 425], [26, 426], [27, 427], [28, 428], [29, 429], [30, 430], [31, 431], [32, 432], [33, 433], [34, 434], [35, 435], [36, 436], [37, 437], [38, 438], [39, 439], [40, 440], [41, 441], [42, 442], [43, 443], [44, 444], [45, 445], [46, 446], [47, 447], [48, 448], [49, 449], [50, 450], [51, 451], [52, 452], [53, 453], [54, 454], [55, 455], [56, 456], [57, 457], [58, 458], [59, 459], [60, 460], [61, 461], [62, 462], [63, 463], [64, 464], [65, 465], [66, 466], [67, 467], [68, 468], [69, 469], [70, 470], [71, 471], [72, 472], [73, 473], [74, 474], [75, 475], [76, 476], [77, 477], [78, 478], [79, 479], [80, 480], [81, 481], [82, 482], [83, 483], [84, 484], [85, 485], [86, 486], [87, 487], [88, 488], [89, 489], [90, 490], [91, 491], [92, 492], [93, 493], [94, 494], [95, 495], [96, 496], [97, 497], [98, 498], [99, 499], [100, 500], [101, 501], [102, 502], [103, 503], [104, 504], [105, 505], [106, 506], [107, 507], [108, 508], [109, 509], [110, 510], [111, 511], [112, 512], [113, 513], [114, 514], [115, 515], [116, 516], [117, 517], [118, 518], [119, 519], [120, 520], [121, 521], [122, 522], [123, 523], [124, 524], [125, 525], [126, 526], [127, 527], [128, 528], [129, 529], [130, 530], [131, 531], [132, 532], [133, 533], [134, 534], [135, 535], [136, 536], [137, 537], [138, 538], [139, 539], [140, 540], [141, 541], [142, 542], [143, 543], [144, 544], [145, 545], [146, 546], [147, 547], [148, 548], [149, 549], [150, 550], [151, 551], [152, 552], [153, 553], [154, 554], [155, 555], [156, 556], [157, 557], [158, 558], [159, 559], [160, 560], [161, 561], [162, 562], [163, 563], [164, 564], [165, 565], [166, 566], [167, 567], [168, 568], [169, 569], [170, 570], [171, 571], [172, 572], [173, 573], [174, 574], [175, 575], [176, 576], [177, 577], [178, 578], [179, 579], [180, 580], [181, 581], [182, 582], [183, 583], [184, 584], [185, 585], [186, 586], [187, 587], [188, 588], [189, 589], [190, 590], [191, 591], [192, 592], [193, 593], [194, 594], [195, 595], [196, 596], [197, 597], [198, 598], [199, 599], [200, 600], [201, 601], [202, 602], [203, 603], [204, 604], [205, 605], [206, 606], [207, 607], [208, 608], [209, 609], [210, 610], [211, 611], [212, 612], [213, 613], [214, 614], [215, 615], [216, 616], [217, 617], [218, 618], [219, 619], [220, 620], [221, 621], [222, 622], [223, 623], [224, 624], [225, 625], [226, 626], [227, 627], [228, 628], [229, 629], [230, 630], [231, 631], [232, 632], [233, 633], [234, 634], [235, 635], [236, 636], [237, 637], [238, 638], [239, 639], [240, 640], [241, 641], [242, 642], [243, 643], [244, 644], [245, 645], [246, 646], [247, 647], [248, 648], [249, 649], [250, 650], [251, 651], [252, 652], [253, 653], [254, 654], [255, 655], [256, 656], [257, 657], [258, 658], [259, 659], [260, 660], [261, 661], [262, 662], [263, 663], [264, 664], [265, 665], [266, 666], [267, 667], [268, 668], [269, 669], [270, 670], [271, 671], [272, 672], [273, 673], [274, 674], [275, 675], [276, 676], [277, 677], [278, 678], [279, 679], [280, 680], [281, 681], [282, 682], [283, 683], [284, 684], [285, 685], [286, 686], [287, 687], [288, 688], [289, 689], [290, 690], [291, 691], [292, 692], [293, 693], [294, 694], [295, 695], [296, 696], [297, 697], [298, 698], [299, 699], [300, 700], [301, 701], [302, 702], [303, 703], [304, 704], [305, 705], [306, 706], [307, 707], [308, 708], [309, 709], [310, 710], [311, 711], [312, 712], [313, 713], [314, 714], [315, 715], [316, 716], [317, 717], [318, 718], [319, 719], [320, 720], [321, 721], [322, 722], [323, 723], [324, 724], [325, 725], [326, 726], [327, 727], [328, 728], [329, 729], [330, 730], [331, 731], [332, 732], [333, 733], [334, 734], [335, 735], [336, 736], [337, 737], [338, 738], [339, 739], [340, 740], [341, 741], [342, 742], [343, 743], [344, 744], [345, 745], [346, 746], [347, 747], [348, 748], [349, 749], [350, 750], [351, 751], [352, 752], [353, 753], [354, 754], [355, 755], [356, 756], [357, 757], [358, 758], [359, 759], [360, 760], [361, 761], [362, 762], [363, 763], [364, 764], [365, 765], [366, 766], [367, 767], [368, 768], [369, 769], [370, 770], [371, 771], [372, 772], [373, 773], [374, 774], [375, 775], [376, 776], [377, 777], [378, 778], [379, 779], [380, 780], [381, 781], [382, 782], [383, 783], [384, 784], [385, 785], [386, 786], [387, 787], [388, 788], [389, 789], [390, 790], [391, 791], [392, 792], [393, 793], [394, 794], [395, 795], [396, 796], [397, 797], [398, 798], [399, 799], [400, 800]];

  return V,EV
end

# //////////////////////////////////////////////////////////////////////////////
function generate_random_bubbles(; n=50)
  store = []
  for k=1:n
    sc = rand()
    scale = S(1,2)(0.25*sc, 0.25*sc)
    transl = T(1,2)(rand(2)...)
    str = STRUCT([ transl, scale, CIRCUMFERENCE(1.)(32) ])
    push!(store, str)
  end
  obj = Plasm.LAR(STRUCT((AA)(STRUCT)(store)))
  V,EV = obj.V, obj.C[:FV]
    return V,EV
  end
  
# //////////////////////////////////////////////////////////////////////////////
function generate_random_lines(; n=400, t=0.4)
  V = zeros(Float64,2,2*n)
  EV = [zeros(Int64,2) for k=1:n]
  for k=1:n
    v1 = rand(Float64,2)
    v2 = rand(Float64,2)
    vm = (v1+v2)/2
    transl = rand(Float64,2)
    V[:,k] = (v1-vm)*t + transl
    V[:,n+k] = (v2-vm)*t + transl
    EV[k] = [k,n+k]
  end
  return V,EV
end
  
# //////////////////////////////////////////////////////////////////////////////
function save_data(V,EV, filename="/tmp/lar.txt")
  open(filename, "w") do f
    write(f, "V = $V\n\n");
    write(f, "EV = $EV\n\n");
    close(f)
  end
end

# //////////////////////////////////////////////////////////////////////////////
function view_data(V,FV,EV)
  
  # disabled since it's another package?
  # using ViewerGL
  # GL = ViewerGL
  # GL.VIEW(GL.GLExplode(V,FV,1.2,1.2,1.2,1));
  # GL.VIEW(GL.GLExplode(V,FV,1.2,1.2,1.2,99,1));
  # GL.VIEW(GL.GLExplode(V,FV,1.,1.,1.,99,1));
end

# //////////////////////////////////////////////////////////////////////////////
function run_test(name, V,EV)
  println("running test ",name)
  save_data(V,EV)
  V,FV,EV = arrange2D(V,EV)
  view_data(V,FV,EV)
  println("done ",name)
end

run_test("generate_debugging_data",generate_debugging_data()...)
run_test("generate_random_lines",generate_random_lines()...)
run_test("generate_random_bubbles",generate_random_bubbles(n=50)...)







  