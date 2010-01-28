function StarMap (elt, size, prop) {
    this.paper = new Raphael(elt, size, size);
    this.size = size;
    var halfsize = Math.floor(size/2);

    this.circle = this.paper.circle(halfsize, halfsize, halfsize);
    this.circle.attr({fill: (prop.circleFill || "#000010")});
    this.constel = this.paper.set();
    this.eqGrid = prop.eqGrid ? this.paper.set() : null;
    this.stars = this.paper.set();
}

StarMap.DEG2RAD = Math.PI/180.0;

StarMap.STARS=[[3,-0.5910663624698816,0.055978128971316],[0,-0.6044947588636527,0.07753736012032166],[2,-0.6180360889141832,0.09288398735791664],[3,-0.682486712751385,0.02718652497282378],[3,-0.5994576209290392,0.00877642375852575],[3,-0.5449220554462875,-0.007168850759277908],[2,-0.7487416157965763,0.12157744609010891],[4,-0.7686340653720335,0.13226047959547785],[3,-0.7458892608738381,-0.04263220645336978],[2,0.018297745742125116,0.25946340040776905],[3,0.0860085228409255,0.27602081076846435],[2,0.15331167524135658,0.3212654017958185],[2,0.2770864371786847,0.3871703329466248],[3,0.12444941437826712,0.3492083497562398],[4,0.10910678718957706,0.3746702253766953],[4,-0.059039078374779376,-0.3425555239467585],[4,0.12855926392780026,-0.2619483862167882],[4,-0.09007505574523632,-0.291774879928412],[3,-0.9714758918790024,-0.4672796554797694],[2,-1.131079570981141,-0.46499073537563945],[3,-1.3113784715888586,-0.5315981300831998],[3,-1.3651691936660344,-0.5662488361088522],[3,-1.1347782968852833,-0.5853319678588527],[3,-1.1636633188534604,-0.5359446467817873],[2,-1.1641478914711458,-0.5264577758277738],[4,-1.7872316022472643,-0.1470396496598434],[3,-1.9676415928844935,-0.12978461960170198],[2,-2.1815750609355646,-0.08206444069389299],[2,-2.5569584588567547,-0.14091114897923918],[3,-2.3548376426731696,-0.22427661974055513],[4,0.29858856558288205,0.07735622418600265],[4,0.3349801803053291,0.07396263228876522],[3,0.23106458124067125,-0.14000297785765742],[2,0.09538005126698286,-0.15826532739345234],[3,0.04241083313841416,-0.07715499673264832],[3,0.15076212576780615,-0.08908865702337936],[3,0.18539361724005288,-0.07152961854195629],[3,0.2480774348936591,-0.09043475029113822],[4,0.32963012690838683,-0.10766759291005974],[4,0.36288892227490294,-0.10397132840554443],[4,0.37418951654278365,-0.1215328744516441],[4,0.34454227012150257,-0.13382216000511213],[6,0.31373193409243155,-0.02598560100511384],[4,0.3626920943577348,0.002866857875600844],[3,0.37215019152432455,0.028248435712653095],[2,0.42005234689495763,0.03570649336344404],[4,0.4134851303231161,0.07788867347566898],[4,0.37622281124281687,0.08849313785449259],[4,0.3538112486168173,0.048849565668346036],[3,0.38883641462601637,0.24248707062139532],[2,0.2847901222251061,0.2076614540224922],[2,0.25545251528274543,0.1836092159387194],[3,0.25287506895802914,0.16998205784868578],[4,-0.8127696379992931,-0.04145634444104267],[5,-0.811653618586904,-0.049827308511687006],[5,-0.7857974857058639,-0.13700983395259317],[4,-0.8800783625091654,-0.12779966842136434],[3,-0.8570762874854032,-0.0720607776689809],[3,2.1457571137581697,-0.31825918745441645],[3,2.1891398532439346,-0.29798383213292806],[4,2.279867761285406,-0.24664217021085424],[3,-2.726201053922747,0.12037876130087984],[-1,-3.2792311174546267,0.1690230679390664],[2,-2.657155608608743,0.24076344760257576],[3,-2.2006362963588137,0.2992069902092346],[3,-2.3855092459236893,0.36783559678853256],[3,-2.9025799568556616,0.3473357402783608],[3,-2.90769373484538,0.27142362128260844],[2,-3.9129806334310056,0.16195201823742908],[4,-4.106905588035603,0.13874195272139309],[5,0.6708167314091483,-0.4137407105230504],[4,0.70204194477507,-0.38249867172486796],[5,0.7069180186410323,-0.3360006829683291],[4,1.9059428225836343,-0.7942964466961625],[4,5.3609371748448105,-0.8186043471936795],[4,-24.9670259619046,-0.8289235611911286],[4,2.2289660565124665,0.2563843620108397],[4,2.1852836388711756,0.18957201172199956],[5,1.9216425553441208,0.2420974917895008],[3,2.2030208058722107,0.15977007370013718],[3,1.885843063269767,0.0803321697748575],[4,2.3918466619101206,0.1038495178810473],[3,-0.5259761589495716,-0.10991399442000088],[3,-0.5177676554292086,-0.1297122891742168],[4,-0.399095445369624,-0.15152778407456233],[4,-0.3584300196810904,-0.14797529079359079],[3,-0.3150852700410198,-0.14643909673584185],[2,-0.2984977089557368,-0.14166732052796485],[3,-0.3475836377589699,-0.1981085229530044],[4,-0.45021516150988594,-0.22417270528452174],[4,-0.4352751868939453,-0.23933245843901582],[1,2.6256833328791287,-0.6965395370220722],[3,4.235992120485922,-0.7007007708295006],[2,5.893359510059905,-0.6296663292197124],[4,5.939021555905657,-0.5839622291216413],[3,8.878462393762092,-0.5654852565776916],[3,6.843546109417638,-0.5640836228594629],[3,4.908528620876219,-0.5627746454439148],[3,4.378669190126615,-0.5929579817509518],[2,2.694142831600548,-0.5689481456693662],[4,2.152068325070008,-0.574572853160198],[1,1.9469845674750108,-0.571657522111942],[-1,1.1103763440087817,-0.4952705540243907],[3,2.5876674378965916,-0.5653919061079064],[3,1.1800006051743364,-0.39588684426512283],[2,1.7637618036409486,-0.3640017549183171],[3,0.25488454356996704,0.6208740909718827],[2,0.1894266253709619,0.5800925658939772],[2,0.12435354241993475,0.5857203978406632],[2,0.08860214233471499,0.5377401673988653],[2,0.020004492797091596,0.5675045914289113],[-1,-2.753426174038431,-0.5871089211713851],[0,-3.61126529548291,-0.5816982635103416],[2,-4.515936399570533,-0.5036734082963681],[2,-3.8827644240805284,-0.4378068718551003],[3,-3.775525889273565,-0.4122063047905832],[3,-4.101509150377092,-0.38861439924455615],[3,-4.105865356034177,-0.38073997808553844],[3,-4.968146911103684,-0.358123517420755],[2,-5.628111688341401,-0.3318038053367029],[2,-3.5254465244036854,-0.3284805726999237],[2,-2.8335845153454104,-0.38544431044694594],[3,-2.426758797511173,-0.38490712736089117],[2,-11.009417025947265,-0.4553035451185675],[3,-16.326357372910326,-0.4687599865045494],[2,-54.828644428060876,-0.4740053815094669],[4,18.14326937854306,-0.5124323039171848],[3,18.909205163445076,-0.6130387968908073],[3,-0.24272383548220292,0.5566071252374268],[3,-0.1546251770496299,0.65190070850187],[3,-0.3427188858317459,0.7075245748483909],[2,-0.3674914262170006,0.6078357297188037],[3,-0.04508238968300348,0.8044800049915815],[4,-6.497809274611108,0.15417483259226064],[4,-6.324288371972348,0.24817928596402405],[4,-16.995487276188175,0.2518193277488242],[4,-13.55604644707336,0.3774374434682392],[2,-8.139838918907207,0.3474370719667896],[2,0.9987818155376496,0.3366611885523757],[1,0.9979490570921412,0.41367623549005134],[0,0.8268614428546872,0.42446484333859885],[3,0.7759520663015089,0.37463971619124914],[2,0.756958233751214,0.2977911751501602],[1,0.8627643392614884,0.25496997128112375],[3,1.1014647426023538,-0.3003594412459458],[4,1.0749623237382093,-0.3166495101984839],[4,0.9893099269851577,-0.3180190091226446],[3,0.961310331553833,-0.3226952746495576],[3,0.9962829126661759,-0.3920478343157894],[2,0.9149233587899878,-0.30643744155603797],[3,0.8816661234884822,-0.31981811624246054],[3,-2.701333483186195,-0.6367583367157773],[4,-2.1038906577425234,-0.5694731243003981],[4,-2.1752060753445854,-0.5634811903323064],[4,7.567579369905259,-0.1610611467709044],[4,9.446576763746693,-0.20186792624049485],[4,13.02702009261164,-0.15555902477191605],[3,11.244133356034157,-0.1296913566786433],[4,12.926154497404664,-0.09505069760298672],[4,19.640651192469445,-0.08575005223400134],[5,115.06117517318258,-0.15079671737296851],[4,30.07067278405561,-0.16152241469905051],[5,-0.82513431406013,-0.34748662529660657],[5,-0.7789176147002204,-0.3379287717386746],[4,-0.7717573081642406,-0.33563480995908895],[4,-0.7453182101526229,-0.3352043933163921],[4,-0.73500549464831,-0.3433931426367533],[4,-0.7331324779775166,-0.3574681644850125],[4,-0.7387838932459008,-0.3688861465856418],[4,-0.7565894525804254,-0.38481546135325745],[5,-0.7802922691056313,-0.3909960132915237],[5,-0.8679276072859017,-0.3610468198861223],[4,-1.9955403163280028,0.2807027751106332],[3,-2.0522255303432853,0.25959998403574996],[2,-1.9765926883393512,0.2374490650468047],[3,-1.8932239953210228,0.23358613446723397],[4,-1.826585183035586,0.2314983777885938],[4,-1.7532901043272249,0.23895328617684247],[4,-1.7195249033814843,0.26655679257969195],[4,-14.267227026456991,-0.1422842693723852],[2,-15.325485523830267,-0.1451275955649905],[2,-28.984744483973344,-0.15428980936191208],[3,-45.261380207385294,-0.19999934537788663],[4,-54.480059604414485,-0.21921280700570125],[2,-13.304540526768015,-0.20705932837573093],[1,-14.68483334924662,-0.5442282949822361],[0,-17.21310304130187,-0.6139904862432507],[1,-9.570205965609803,-0.5737341571044958],[2,-30.2522904688559,-0.562883691494912],[3,-0.709672536808893,0.5025995048801761],[3,-0.6691077992441574,0.48481351809336043],[2,-0.6219185537574116,0.41555039015833456],[2,-0.514403528218701,0.3665090236148914],[1,-0.46250997227335416,0.41708258020984873],[2,-0.4499235481776471,0.30543921522643386],[3,-0.38152193140231283,0.27007426566559556],[4,-0.3053967868922274,0.25622828525456376],[3,-0.5881441579697033,0.3161000151179943],[3,-0.6659044894470283,0.24895442255447198],[4,-0.4844615793169225,0.09896134654469271],[3,-0.4728369213548544,0.12806040425095658],[3,-0.46727052989688006,0.1397582813768529],[4,-0.4487397336072945,0.14165135915461222],[4,-0.45715206383424856,0.13231557862682763],[4,0.9356756575232484,-0.6461313712973283],[4,0.9745750370537413,-0.6138939273910802],[3,0.8910664418444603,-0.6066938175756883],[3,0.6808607250522875,-0.5210663455011186],[4,0.6249448151445668,-0.4822039210541508],[3,-1.0286499385688204,0.5415153422727457],[2,-1.0149183587771846,0.48222388487530427],[2,-1.1381110577755844,0.49099223002910714],[4,-1.1289801247445472,0.5224852033595456],[3,-0.7247021835008325,0.6702167782111405],[3,-0.612290190818346,0.70369652391748],[3,-0.9121701444517182,0.7363708709552264],[3,-1.2527504168606949,0.6458728073864836],[2,-1.5399255249482215,0.5951038149193899],[4,-1.7156305690913738,0.5607658001591382],[3,-2.0856453461489597,0.5653814038751053],[3,-3.593991808540786,0.6294388512682881],[3,-13.664713843449398,0.6974566973650491],[3,16.00866817182616,0.6915439393734262],[4,-1.6200953775168891,-0.5164746343431237],[4,-1.5708742186653226,-0.4679595199793331],[4,-1.5166937043508946,-0.4405808309324826],[4,-1.7043361593645614,-0.4581495061750619],[0,0.21646462314573323,-0.5456328053470814],[3,0.2584947086974413,-0.4835221452651111],[3,0.30694376891035563,-0.4824730549850775],[4,0.33213502535254213,-0.44213530655372313],[4,0.36347287742167494,-0.3928184427791113],[4,0.3656154194783631,-0.36254094758894384],[2,0.40977803972269544,-0.3669847558309095],[4,0.4660077658124897,-0.39463190897466743],[4,0.5446434101412764,-0.340623804819304],[4,0.5470720035001884,-0.32685146444141294],[3,0.6306286587930789,-0.3038073437643784],[3,0.6495195638561677,-0.30589278529608116],[3,0.6858403236404361,-0.2732155514066401],[4,0.5397110734883597,-0.20571107294169097],[4,0.503460149580211,-0.19105701616081996],[3,0.4650231888411586,-0.1921890451556208],[4,0.42034388031162895,-0.2091325383811606],[4,0.3766216560443497,-0.16351122052698383],[3,0.4051095820834275,-0.07780270955708972],[4,0.455288398420319,-0.07711997211186142],[3,0.5011457603334452,-0.08272724675590652],[3,0.5296103734056369,-0.08542432012816945],[3,0.6883061523728624,-0.029263978263997865],[4,0.7182499125468553,-0.028409597117152355],[4,0.7429861209317123,-0.04762021898821316],[2,0.7949037476534421,-0.04441518571284389],[4,0.7995275506205445,-0.07654268155550127],[3,0.6943103714384528,-0.1254746736473147],[4,-0.6338575908192213,0.15370142806670067],[3,-0.6146398638992528,0.16316726234633006],[4,-0.6367745222534738,0.15850912922452232],[3,-0.580973274734083,0.17176037692730153],[5,-0.5624425768364378,0.17624488993304574],[4,0.3865932448925073,-0.2905863153284031],[3,0.4454151779671595,-0.2585161586512913],[1,1.1797421504480354,0.1440969768488954],[4,1.3277010988155271,0.1814630668098721],[3,1.428965057731631,0.19422038172177764],[3,1.4155794346359925,0.14535376552412232],[3,1.2200797103055612,0.1130168416338202],[4,1.539297876091095,0.2391191211879875],[3,1.6038189117707136,0.21619056752703772],[1,1.6106584285606418,0.24957189513993086],[3,1.4668245545156269,0.24745910238941335],[4,1.370964445386506,0.27024467130136887],[1,1.5296168223674769,0.2856939320254948],[3,1.2616238355733296,0.30536208750643773],[3,1.21273137261063,0.22289577160921925],[4,1.1350917084576044,0.17823668949315613],[2,1.1055568976429349,0.19903802498775858],[3,1.0671210139678633,0.1989742707596507],[4,1.0181415861933913,0.2058492868007382],[4,0.4906226011749689,0.576656224964005],[4,0.5698809365398672,0.6136684194388505],[4,0.7469059071856639,0.653667226164062],[4,0.5496350975175164,0.7176783136302085],[5,0.8486730375777197,0.827724174401644],[4,1.2695837229283415,-0.10544474862357832],[4,1.3255893834567989,-0.13727866045261217],[4,1.2807312385273661,-0.1499350005584604],[-2,1.2193397409032765,-0.146892732784633],[3,1.3211857084734298,-0.21103659133173497],[1,1.3538388793050586,-0.23448534770980559],[4,1.3942838600057272,-0.23798173565970016],[2,1.4556559790229116,-0.26144166508922545],[1,1.295166646652731,-0.25835774888183327],[3,1.313396147042308,-0.24872413120842232],[6,1.2697324130059642,-0.21190312293313193],[3,1.1743901365764335,-0.16963743768884473],[4,1.1659802530891366,-0.2031327189609117],[1,1.1042970361878692,-0.1579901156309805],[3,1.2450897569967538,-0.29155379163057804],[3,1.0928063061356137,-0.26854204940194654],[1,-4.183694534127367,0.4590323022002961],[2,-5.40053958276852,0.5197400452405324],[1,-8.444657356909657,0.5312600170972344],[3,-29.706282949099347,0.5433240997733226],[1,8.104962543894972,0.5979084490649779],[2,7.838762426595809,0.5359969535117323],[2,74.27335024533475,0.5061744438337907],[3,32.85296475628549,0.44292329434154293],[3,9.069549781511727,0.40909606034469137],[3,4.379343557791358,0.393048496211488],[3,4.621770012655005,0.37886044967842425],[4,3.490460952879728,0.5102322754246541],[3,3.007689980883454,0.484265361318603],[3,2.4692835135241302,0.4364388973854095],[3,2.4025526156145585,0.44567158353466396],[3,3.4587714682879476,0.5662233650989948],[3,2.0307826396841024,0.5857400774365584],[3,2.9784492588941194,0.6135427630593798],[4,-0.11641245893072423,-0.3991651564481573],[3,-0.20057005551518653,-0.3989146029150134],[1,-0.24879193594385876,-0.43440380533354794],[2,-0.1703369898421291,-0.4336151197603359],[3,-0.10872932944154157,-0.4167374929645158],[4,-0.12969719967781587,-0.4959052119138402],[3,-0.15714720730549356,-0.48037255679460145],[4,-0.2537003071143163,-0.35946083875677387],[3,-0.2821997077375407,-0.3381391956741466],[3,-1.0938698066159893,0.4245399270773126],[3,-1.016484984099258,0.33702807054986555],[4,-1.1724911405170424,0.3360119102208968],[4,-1.20412415571969,0.33742483369026505],[3,-1.2182453288120547,0.3327444431362604],[3,-1.4091173754903203,0.3533572670665723],[4,-1.467986584828158,0.3882449144207571],[3,-1.5716465408141462,0.4277049263677354],[4,-1.7978458933062855,0.3883763893819727],[2,-1.4196192821843037,0.2829892522723876],[2,-1.4950447155191295,0.18976235199578817],[3,-1.5552650194931272,0.16871568887424906],[3,-1.3015190876964142,0.276632418862551],[4,-1.1365834084916755,0.23188507187064258],[3,-1.2183295263416904,0.22022646670413082],[3,-1.0608866300634139,0.24675736232639356],[3,-1.0098053211601583,0.2609271444935179],[3,-0.9676197182373801,0.2564072868494057],[3,0.6188228582196991,-0.3868093609238853],[5,0.36560090458956684,-0.5155852097130601],[5,0.4234840441818204,-0.5743013962289826],[4,2.184501422533422,0.02966764910902066],[4,2.12941761933976,0.02916813683392429],[4,2.116195188707679,0.04981618033526751],[3,2.230010845551579,0.056074186311002645],[4,2.2517594432295382,0.05098927689483171],[3,2.347230240084611,0.05193112177815114],[3,2.6457809165745205,0.020205145790027416],[4,2.9882999778049224,-0.01033842782385393],[4,2.9279990558789963,-0.02416840513309485],[1,2.895743198812154,-0.07570540449867935],[4,3.472485225894189,-0.13029042253883377],[3,4.109565290258336,-0.1082271972658652],[3,4.812516454276552,-0.14798903241179667],[3,6.46183154746511,-0.14226883062790213],[3,16.959605877949866,-0.285400300971353],[4,64.63864478489732,-0.30485423810829904],[4,-6.587437422800686,-0.20452464652127159],[2,-5.750413702047388,-0.20501053116842832],[2,0.05602547658767529,-0.7990761203162127],[3,0.5408028675727768,-0.7568325433734726],[4,0.3629503516569132,-0.6779032469051871],[4,0.3195022400641515,-0.6829138792124844],[2,0.265062385680654,-0.5957639833550494],[4,-0.3643044528522256,-0.5034863794879814],[3,-0.4727872035581259,-0.4378421850577145],[3,-0.42754284732704967,-0.5595007851887516],[4,-0.23093824538815194,0.34187677697935553],[4,-0.19780652918004252,0.395154861781296],[4,-0.19997758287102693,0.4421673186997568],[4,-0.2113788296183679,0.4607566513404474],[4,-0.2135588072253118,0.4902141774684407],[3,-0.19598938730345325,0.4693119175990238],[3,1.0670152907288362,-0.05481212119556019],[3,1.1343667039603063,-0.06145212363978122],[4,1.3755257597307202,-0.004300348079844537],[4,1.1094840430992154,0.04010155415436241],[5,1.0399031810346981,0.02181769209310324],[4,1.8095768782446309,-0.02604424809333778],[3,1.5791568667026943,-0.08354247402104747],[4,1.027225971743473,-0.13107802284379072],[3,0.9844369106061583,-0.12427334779670451],[3,0.9446444677731783,-0.130072086162003],[2,0.8875693322986085,-0.15679513056644698],[3,0.813153425040692,-0.14236941649653925],[3,0.9628225623956501,-0.18423365746295203],[3,0.9344186523489613,-0.19843601101689926],[2,0.8702273355959007,-0.18316661304413578],[3,0.7864310751170243,-0.1977410537380101],[4,0.8375224545413993,-0.1154985797604645],[4,0.8142419687293065,-0.11341663171262911],[4,0.8108615526667441,-0.10394980994105997],[5,0.8390396516701312,-0.10788963885821831],[2,41.908077269659955,0.1278575950218267],[3,9.983588981660002,0.13546979090946173],[1,4.024779976195682,0.10481481151776381],[3,3.986046444367985,0.14733425537296135],[2,4.509281709532513,0.1749043853742648],[2,9.954294743932632,0.1810460605986467],[3,4.361422092351806,0.20724775789588992],[3,3.509539324574528,0.23093331438830342],[2,3.3187449569364524,0.21049864643281851],[3,-1.8137496648974876,-0.30217593055343905],[5,-1.6703023730158335,-0.3322248449453553],[3,-1.7309837061434057,-0.3482039682275471],[2,-1.9717369570911143,-0.3755436134252915],[3,-2.1278516181390326,-0.37038215736030955],[3,-2.1226163467166175,-0.32744126848331445],[2,-2.4362449464339275,-0.3952605681402564],[4,-1.9410440934756836,-0.3895568737572464],[3,-2.2422605931043353,-0.4888033339495044],[2,-2.711902128313833,-0.4388458984106622],[4,-2.787381092502394,-0.4602209675545372],[4,-3.0285949115523674,-0.4180964063315676],[3,2.7673280561382154,0.30948023676071135],[3,2.7260747427712864,0.33268381631525595],[4,2.514961605364571,0.34874819070247126],[3,2.423868230758549,0.3816973303432654],[4,1.9503455946292156,0.39581065682102773],[4,1.47362826430274,0.4579586446780175],[4,1.2873158751621747,0.5591452472999854],[4,1.0895088972484546,0.565898411255877],[0,-0.8505538005062213,0.35198878501991054],[4,-0.8215108177178154,0.34047692492050374],[3,-0.8022980240583429,0.29965940391067547],[3,-0.7709947952666509,0.2932689632906395],[4,-0.7865525355492848,0.3336117102620954],[4,4.869062014844691,-0.2779618358600169],[5,3.6956576062593127,-0.3238677574565626],[4,-0.3690744014538285,-0.28837507475460195],[4,-0.4109176870101295,-0.28918135273981316],[4,-0.4400938024200275,-0.3036294822403994],[3,-9.870261781759092,-0.6758818579519346],[3,31.84143694224488,-0.6584879362603411],[3,-14.09361888875574,-0.7283169879012555],[2,-12.299665451254102,-0.689023840048168],[3,-0.3117654853196658,-0.8010001466489752],[4,-0.16272146794043382,-0.859856628470779],[4,-3.012055327601529,-0.8951692243339445],[3,-2.6074337904966547,-0.8249919808328766],[3,-1.4724378775178075,-0.8228272550336443],[4,-1.4078501713967713,-0.8028226244435085],[2,-1.115822718826233,0.1100528049482278],[2,-1.0748442679813022,0.0398749814054531],[2,-1.2438732227898615,-0.13809574037422978],[3,-1.3170264104385763,0.08199582176812335],[3,-1.5824498433782113,-0.040973634788280165],[2,-1.447157683735577,-0.09247808164498404],[4,-1.1332009369244047,-0.21221527210180327],[1,0.9193810835108198,-0.016953819300707627],[1,0.9012479131131417,-0.010489091807840464],[2,0.8847511010182809,-0.0026100882246509628],[5,1.0151943754043202,0.17354442951332844],[4,1.0535036755875529,0.1246346995329575],[4,1.0335972669093179,0.1295980126749461],[4,0.9757941391648747,0.17881521903290135],[4,1.0104530044462814,0.0843886336883426],[0,0.979150500162884,0.06472875862238704],[2,0.9479553088240286,-0.08458404841270446],[0,0.8189928374642708,-0.07169527022558621],[1,0.858357565443754,0.055468667735297],[3,0.8970044026913446,0.08690976083106343],[3,0.7326790628998563,0.060823184445203585],[3,0.7372840628258088,0.04895278671826183],[5,0.7446281220705693,0.02189223848897396],[4,0.7623074888880964,0.014958892689280068],[4,0.7352854378856701,0.07782591063625692],[4,0.7497815252026184,0.08881790459228918],[1,-0.5050063697994688,-0.5399638397862041],[4,-0.3481352874255399,-0.6415960725164714],[3,-0.45320139026456147,-0.6519320616515419],[3,-0.5523211058204007,-0.6516347098774601],[3,-0.575635634254005,-0.7387488695606388],[4,-0.827876594535581,-0.7189398332981486],[4,-0.777950794721592,-0.664820302039397],[4,-0.7946623727641206,-0.6030905503736804],[4,-0.9034620616870325,-0.5948654849515161],[4,-0.9632471838261318,-0.6208487313063261],[3,-1.0642709934937096,-0.6336860686208741],[2,0.028884138466089625,0.1332829592919069],[2,-0.12110133916179536,0.13347616460740394],[2,-0.1232889999938922,0.2500924246666659],[2,-0.1695811746825864,0.2700201960116974],[5,-0.24648440253343953,0.29785162036618934],[3,-0.153910865973201,0.21805066489319644],[3,-0.16167197430853084,0.20859832679797272],[3,-0.25162814880749773,0.22485586697223736],[4,-0.3041926301416696,0.22760745681007777],[4,-0.16131402291832495,0.10664043366665808],[3,-0.1730415996215229,0.09480423010554768],[3,-0.24424283581044162,0.05413871366865765],[2,-0.3052878801628682,0.08638965620591903],[3,1.2359570200455734,-0.6001717115249617],[4,0.9565657548036989,-0.5335744929075613],[3,0.9460031686950645,-0.477690345849548],[3,0.5326018007533848,0.28946986616007825],[2,0.5604052034906508,0.28564694394635226],[3,0.5743434803362238,0.3229047278994097],[2,0.5711236712104304,0.36407176415307857],[3,0.5287040012803863,0.44301014891004337],[1,0.47783944067164635,0.46483434373757704],[2,0.426527170858933,0.5041120708928263],[3,0.3906291971045579,0.5305405815500475],[2,0.43525040824977523,0.37344361327767983],[3,0.42750331333277863,0.35255324218851186],[4,0.39034096165934484,0.34744357755799354],[4,-0.3880219642876439,0.08864906391291828],[4,-0.3776693340362825,0.08755653286535786],[5,-0.35683859987095423,0.05950826967192577],[3,-0.3743243036038191,0.04583013390230107],[0,1.5644754959345601,0.04565028176514712],[2,1.4766551170141213,0.07246518094478493],[3,6.824590671349204,0.3077888559863139],[4,4.908845814601061,0.3317592180392886],[4,3.9895929283032383,0.3176479813083697],[4,3.037562854451136,0.3287601975032615],[4,-0.6722670603824541,0.2186310400966725],[4,-0.5741542383897628,0.24704498468323527],[1,0.34380202376024543,0.9872385236003124],[4,-1.129226766157169,0.9421270552340677],[4,-1.3892832174762781,0.8698577045995147],[4,-1.8801549562187612,0.8068190447683354],[4,-1.5886772492999603,0.7778441692162661],[3,-2.135630785964772,0.7243312573418024],[2,-2.559806840472796,0.7556834796678599],[3,0.1503083380950195,-0.5232967131771157],[3,0.1451826847620237,-0.4318944050076719],[3,0.05722610408651366,-0.40078057459239463],[3,0.20174819702885813,-0.45649666534518124],[4,0.25314800356304173,-0.42759318005558783],[3,0.1952083994892264,-0.39711692548265914],[2,0.057399430479194755,-0.3869222208545399],[3,0.02053013175431177,-0.4218714993718347],[5,0.13791237449164173,0.28489840583054465],[4,0.16229810208820414,0.21788685512430783],[4,0.17512667271435775,0.24251755392681762],[3,0.202278482559029,0.1347239808327503],[4,0.23407181614754505,0.08008606101929966],[3,0.27274070265622713,0.02412303847898802],[4,0.2529358591548045,0.027823173747381132],[4,0.22497414851763686,0.04792502445103683],[4,0.19932603343180932,0.05366736519318388],[4,0.1381928355530373,0.06896288750965829],[5,0.10574190809031242,0.063790129090714],[5,0.04496756700654738,0.07159535763048304],[4,-0.0015051448833216087,0.059967885182598185],[4,-0.04377671999961173,0.0491474387307661],[4,-0.03918523362816187,0.015538321018449004],[4,-0.07226914682908568,0.010959689146337614],[3,-0.09373940737510546,0.02865082954973421],[4,-0.06999379628768006,0.055725704531988295],[1,-0.13687828942013697,-0.26441529162072236],[4,-0.1748531129449383,-0.24048135984871546],[5,-0.2659932775382187,-0.2535379289400644],[5,-0.29683963008821734,-0.2763690682528],[4,-0.24437161683598033,-0.2919321086216597],[4,-0.1955015665040885,-0.2900161697126917],[4,-0.1406559721559613,-0.291849377746974],[3,1.5835254902662053,-0.7346554184706467],[3,1.3560391796861069,-0.7067195603987579],[4,1.803404535491563,-0.6823728807035353],[3,1.407335174606669,-0.6739649349850879],[3,1.9811391404473222,-0.6511027067453757],[4,2.4511382132041133,-0.6543298394515173],[2,1.7998323937767133,-0.21533488534990525],[5,1.6400991617702894,-0.2208913262497695],[2,1.409367307604384,-0.33554157695915127],[2,1.2456156668393377,-0.47285158503299063],[3,1.4911860299002737,-0.3969569207657507],[3,0.6200983201483689,-0.606504262309078],[4,0.6263407114549323,-0.5692541168017488],[4,0.5737071526460871,-0.5937583488944347],[3,0.5322485075778475,-0.6347060786788046],[2,-0.9123509526202049,-0.2663414545906701],[2,-0.8848373169927319,-0.22555457729256703],[3,-0.9258935915725862,-0.33228049081280187],[1,-0.8997228447021574,-0.30940172833407403],[2,-0.9749707079141149,-0.27191562336283404],[4,-1.0558062824037244,-0.2477598864719286],[2,-0.7583124864210615,-0.26682819143147957],[3,-0.8182868146896022,-0.23999366018715199],[3,-0.9416792326915956,-0.18587012219313875],[3,-0.7435438295956776,-0.2462696693661238],[2,-0.7838702584944932,-0.23359605738321085],[3,-0.7752266232714913,-0.18630164875447683],[3,-0.7512229749095851,-0.19203894376347666],[4,-0.7079254842061009,-0.16691995208346957],[3,-0.694782055315807,-0.15701831595021282],[5,-0.6493224225855657,-0.2191235140850002],[4,-0.5696449869763214,-0.246640531393633],[4,-0.5781178475996069,-0.3179511476414547],[4,-0.591217723211698,-0.38254707601027715],[3,-0.6876476888778357,-0.37006645080600514],[4,-0.6897969329848236,-0.4121666535813489],[1,-1.1223321116741773,-0.3356023601133513],[2,-1.0794857358212455,-0.35441247332084924],[2,-1.055694223490702,-0.3652256874683039],[1,-1.104207269750305,-0.3938885216385888],[3,-1.2340547047211654,-0.39631661696200554],[4,-1.33918974696777,-0.3874931641017439],[3,-1.3522200871904964,-0.344789661339411],[2,-1.362762651316903,-0.30852506066829216],[2,-1.4558072538387425,-0.2513307245923666],[1,-1.5007998935161815,-0.23484202664625006],[2,-1.7291435532515582,-0.20001592298017062],[2,-1.7421130425926057,-0.231916692983512],[2,-1.6855551314826684,-0.17457644295841876],[3,-1.826337546229532,-0.029942573546032553],[3,-1.8150999106935741,0.03909415054462012],[2,-1.8780906415451837,0.05613208301674512],[3,-1.9753446155726269,0.09222907618158759],[3,-1.8592721399102174,0.13540014203120088],[3,-1.7634531658733315,0.13755850107730613],[4,-1.8346996270298244,0.15965298344766657],[4,-0.7805152285330654,0.03669917727009839],[3,-0.9110661027851723,-0.025287547081553362],[4,-1.084566401769194,-0.11283228502946456],[3,-1.1029090469007963,-0.13519115546451482],[4,-1.1873946187640463,-0.1125822549960283],[5,5.0441135263967585,-0.005558687327971446],[4,4.008480398327549,-0.0032430761934696952],[5,0.8842496422512588,-0.7860984393765768],[5,0.7102223568549931,-0.7123864581044935],[4,0.7075320913827426,0.20306133285212533],[3,0.6638185967922692,0.16896249486011058],[0,0.6870240841292656,0.14508008331918584],[2,0.9069225579959002,0.18662644583249705],[3,0.6364332287466583,0.13722923621951286],[3,0.6461054259010742,0.15429515331448598],[3,0.5793310088177449,0.1094328153767302],[3,0.47915720658822053,0.07895692487917032],[3,0.663962190170179,0.13939257384011233],[4,0.6540346414657421,0.15774033156857578],[3,0.5462432000850487,0.21304480491355895],[4,-0.881506219885216,-0.45646453466740805],[3,-0.8887242742357389,-0.42414869633079194],[2,-0.22512386564370243,-0.5803735859905252],[3,-0.09314081990797124,-0.5570036045120967],[4,0.043730675913326654,-0.6355699404794433],[4,0.06892428280069794,-0.6122977857583141],[5,0.30558337159937304,0.2996241769986332],[3,0.29038970773017214,0.3151778250172617],[3,0.25183600476930396,0.2640190721925466],[1,-1.372194821964246,-0.6876361054207837],[2,-2.157870892923492,-0.683170370550704],[2,-1.7751956523294765,-0.6179712667845698],[2,-0.33566611450540607,-0.04865584899866612],[2,-0.2544697690077704,-0.0027910206926115063],[3,-0.2179095097839653,-0.012107531897251741],[3,-0.20156735920964514,-0.00017508306612997157],[4,-0.18679382442893974,-0.0010241702874522129],[3,-0.1480808957572566,-0.06624181463181524],[4,-0.09653584509810018,-0.07947177636197722],[3,-0.08095941427036027,-0.17723002699712198],[4,-0.22895526353404586,-0.06802619578031177],[4,-0.25295364080625804,-0.12162908470399908],[4,-0.19744541302533028,-0.09345277268943221],[4,-0.15482650024841169,-0.11917673783689275],[3,-0.14354403526345436,-0.13894612876412427],[3,-0.11074096157521594,-0.18689640028984653],[3,-0.4460770509775242,-0.0830557249131587],[4,32.40480639298109,0.057045131576573226],[5,-24.535255808605413,-0.006869242077650569],[2,-10.970468739569897,-0.012650142603547517],[0,-5.318204923146477,-0.0977094162519204],[4,-3.3518773080987394,-0.0898989629888957],[4,-3.270475537572082,-0.052403633526946415],[3,-2.691445166294534,-0.049410535299839825],[3,-4.771341471619934,-0.005200588627928909],[4,-3.679140206568248,0.013479848041699074],[3,-2.6351187642155445,0.016520505955691808],[3,-8.202321255779847,0.029658330906551154],[2,-7.326405100765288,0.09592884641885469],[1,1.818359703325271,-0.43830971113230316],[3,2.148096265160659,-0.49773172992425685],[1,2.2032520927904327,-0.517338890258013],[2,2.787407366616979,-0.5206856913868292],[3,3.6324152172724213,-0.5157826722194386],[2,6.205807477745477,-0.4601614038630822],[3,5.482535726646758,-0.44758989303846586],[3,4.27770503116962,-0.38508668651864225],[3,2.9607926112211636,-0.3685912451543431],[2,2.538587596391352,-0.3982781567749169]];

StarMap.CONSTELLATIONS=[[[3,8],[6,7],[3,6],[5,4],[3,4],[1,3],[1,2],[0,1],],[[13,14],[11,13],[12,11],[10,11],[9,10],],[[17,15],[16,17],[15,16],],[[24,18],[23,24],[22,23],[21,22],[20,21],[19,20],[18,19],],[[29,26],[28,29],[27,28],[26,27],[25,26],],[[48,44],[31,48],[47,31],[46,47],[45,46],[44,45],[43,44],[42,43],[42,39],[41,32],[40,41],[39,40],[38,39],[37,38],[36,37],[35,36],[33,35],[33,34],[32,33],[30,31],],[[51,52],[50,51],[49,50],],[[57,53],[56,57],[55,56],[54,55],[53,54],],[[59,60],[58,59],],[[68,69],[62,68],[67,62],[66,67],[65,66],[64,65],[63,64],[62,63],[61,62],],[[71,72],[70,71],],[[74,75],[73,74],],[[79,81],[79,80],[77,79],[77,78],[76,77],],[[84,90],[83,89],[88,84],[85,88],[86,87],[85,86],[84,85],[83,84],[82,83],],[[101,105],[102,104],[103,100],[103,99],[101,102],[100,101],[98,99],[97,98],[96,97],[95,96],[94,95],[93,94],[92,93],[91,92],],[[109,110],[108,109],[107,108],[106,107],],[[126,127],[125,126],[124,125],[123,124],[114,123],[121,122],[116,121],[117,120],[118,119],[117,118],[116,117],[115,116],[114,115],[113,114],[112,113],[111,112],],[[132,130],[129,132],[131,128],[130,131],[129,130],[128,129],],[[134,135],[133,134],],[[136,137],],[[143,138],[143,142],[141,142],[140,141],[139,140],[138,139],],[[149,150],[147,149],[147,148],[146,147],[145,146],[144,145],],[[151,153],[151,152],],[[161,156],[160,161],[159,160],[158,159],[157,158],[157,154],[156,157],[155,156],[154,155],],[[162,171],[169,170],[168,169],[167,168],[166,167],[165,166],[164,165],[163,164],[162,163],],[[177,178],[176,177],[175,176],[174,175],[173,174],[172,173],],[[184,180],[182,184],[182,183],[181,182],[180,181],[179,180],],[[187,188],[185,186],],[[197,198],[192,197],[195,196],[194,195],[192,194],[192,193],[191,192],[190,191],[189,190],],[[203,200],[202,203],[201,202],[200,201],[199,200],],[[207,208],[206,207],[206,204],[205,206],[204,205],],[[221,222],[220,221],[219,220],[218,219],[217,218],[216,217],[215,216],[214,215],[213,214],[209,213],[212,209],[211,212],[210,211],[209,210],],[[226,223],[226,224],[225,226],[224,225],[223,224],],[[253,254],[252,253],[251,252],[250,251],[249,250],[248,249],[247,248],[246,247],[245,246],[244,245],[243,244],[242,243],[241,242],[240,241],[239,240],[238,239],[237,238],[236,237],[235,236],[234,235],[233,234],[232,233],[231,232],[230,231],[229,230],[228,229],[227,228],],[[258,259],[256,258],[256,257],[255,256],],[[260,261],],[[277,278],[276,277],[274,276],[274,275],[271,274],[271,273],[271,272],[270,271],[267,270],[267,269],[267,268],[264,267],[265,266],[264,265],[263,264],[262,263],],[[282,283],[282,281],[279,282],[280,281],[279,280],],[[286,284],[299,292],[292,298],[295,287],[295,297],[295,296],[294,295],[293,294],[293,289],[292,293],[290,291],[289,290],[288,289],[287,288],[286,287],[285,286],[284,285],],[[317,304],[316,317],[315,316],[311,315],[312,314],[312,313],[311,312],[305,311],[308,310],[308,309],[307,308],[306,307],[306,303],[305,306],[304,305],[303,304],[302,303],[301,302],[300,301],],[[325,326],[320,325],[321,324],[321,323],[322,318],[321,322],[320,321],[319,320],[318,319],],[[339,331],[343,342],[343,344],[342,343],[340,341],[339,340],[336,339],[337,338],[336,337],[332,336],[334,335],[333,334],[332,333],[331,332],[330,331],[329,330],[328,329],[327,328],],[[346,347],[345,346],],[[364,365],[363,364],[362,363],[361,362],[360,361],[359,360],[358,359],[357,358],[356,357],[355,356],[354,355],[353,354],[352,353],[352,348],[351,352],[350,351],[349,350],[348,349],],[[369,370],[368,369],[367,368],[366,367],],[[373,371],[372,373],[371,372],],[[379,376],[378,379],[377,378],[376,377],[375,376],[374,375],],[[385,386],[382,385],[383,384],[382,383],[381,382],[380,381],],[[396,399],[397,398],[395,391],[391,397],[391,396],[390,394],[394,395],[393,394],[392,393],[390,392],[390,391],[389,390],[388,389],[387,388],],[[405,401],[407,408],[406,407],[404,406],[405,400],[404,405],[403,404],[402,403],[401,402],[400,401],],[[418,415],[418,420],[417,419],[417,418],[416,417],[412,416],[413,415],[413,414],[412,413],[411,412],[411,409],[410,411],[409,410],],[[427,428],[426,427],[425,426],[424,425],[423,424],[422,423],[421,422],],[[433,430],[432,433],[431,432],[430,431],[429,430],],[[434,435],],[[437,438],[436,437],],[[442,439],[441,442],[440,441],[439,440],],[[445,443],[444,445],[443,444],],[[447,448],[446,447],],[[451,455],[454,451],[453,454],[452,453],[449,452],[451,450],[449,450],],[[461,463],[473,474],[469,473],[471,472],[470,471],[469,470],[467,469],[468,464],[467,468],[458,467],[466,458],[465,466],[456,465],[464,456],[463,464],[460,463],[461,462],[460,461],[459,460],[457,458],[456,457],],[[484,485],[484,482],[483,484],[482,483],[481,482],[481,478],[480,481],[479,480],[478,479],[478,475],[477,478],[476,477],[475,476],],[[488,487],[9,486],[9,488],[497,498],[496,497],[495,496],[487,495],[493,494],[492,493],[491,492],[488,491],[489,490],[488,489],[486,487],],[[500,501],[499,500],],[[511,512],[510,511],[507,510],[508,509],[507,508],[506,507],[505,506],[504,505],[503,504],[502,503],],[[516,513],[515,516],[514,515],[513,514],],[[517,518],],[[521,519],[521,522],[520,521],[519,520],],[[523,524],],[[531,528],[530,531],[529,530],[528,529],[527,528],[526,527],[525,526],],[[539,534],[538,539],[534,538],[537,534],[533,537],[536,533],[535,536],[533,535],[534,532],[533,534],[532,533],],[[557,553],[556,557],[555,556],[554,555],[553,554],[552,553],[551,552],[550,551],[549,550],[548,549],[547,548],[546,547],[545,546],[544,545],[543,544],[541,543],[542,541],[540,542],[540,541],],[[563,564],[562,563],[561,562],[560,561],[559,560],[558,559],],[[570,567],[569,570],[567,569],[567,568],[567,565],[566,567],[565,566],],[[105,571],[575,105],[574,575],[104,574],[573,104],[572,573],[571,572],],[[579,576],[578,579],[577,578],[576,577],],[[598,600],[598,599],[597,598],[596,597],[595,596],[589,595],[593,594],[592,593],[591,592],[590,591],[590,587],[589,590],[586,589],[581,588],[587,581],[587,580],[586,587],[583,586],[580,583],[584,580],[584,585],[583,584],[582,583],[580,581],],[[610,613],[610,612],[610,611],[609,610],[608,609],[607,608],[606,607],[605,606],[604,605],[603,604],[602,603],[601,602],],[[624,625],[623,624],[622,623],[621,622],[620,618],[619,620],[618,619],[617,618],[616,617],[615,616],[614,615],],[[626,627],],[[628,629],],[[635,640],[639,635],[631,639],[638,634],[632,638],[632,631],[636,637],[634,636],[634,635],[632,633],[630,631],[143,630],],[[641,642],],[[644,646],[644,645],[643,644],],[[649,647],[648,649],[647,648],],[[652,650],[651,652],[650,651],],[[667,653],[665,666],[664,665],[663,664],[661,663],[661,662],[654,661],[659,660],[658,659],[657,658],[656,657],[655,656],[654,655],[653,654],],[[678,670],[678,679],[675,678],[676,677],[675,676],[671,675],[673,674],[672,673],[671,672],[670,671],[669,670],[668,669],],[[689,680],[688,689],[687,688],[686,687],[685,686],[684,685],[683,684],[682,683],[681,682],[680,681]]];


function stereographicProjectPoints(arr, lam1, phi1, rad) {
    function sinSum(cosa, sina, cosb, sinb) {
        return cosa*sinb+sina*cosb;
    }

    function cosSum(cosa, sina, cosb, sinb) {
        return cosa*cosb-sina*sinb;
    }
    var DEG2RAD = StarMap.DEG2RAD;
    lam1 *= DEG2RAD;
    phi1 *= DEG2RAD;
    var len = arr.length, i;
    var res = Array(len);
    var cphi = Math.cos(phi1), sphi = Math.sin(phi1);
    var clam = Math.cos(lam1), slam = Math.sin(lam1);
    for (i = 0; i < len; ++i) {
        var star = arr[i];
        var mag = star[0], re = star[2], de = -star[1];
        var t2c = re*re, t2l = de*de, t2c1=1+t2c, t2l1=1+t2l;
        var cosc = (1-t2c)/t2c1, sinc = 2*re/t2c1;
        var cosl1 = (1-t2l)/t2l1, sinl1 = 2*de/t2l1;
        var cosl = cosSum(cosl1, sinl1, clam, slam), sinl = sinSum(cosl1, sinl1, clam, slam);
        var k = rad / (1.0 + sphi * sinc + cphi * cosc * cosl);
        var x = k * cosc * sinl, y = k * (cphi * sinc - sphi * cosc * cosl);
        res[i] = [mag,
                  x,
                  y,
                  x*x + y*y < rad*rad];
    }
    return res;
}

StarMap.prototype = {
    setPos: function (lat, lon) {
        var ortho = stereographicProjectPoints(StarMap.STARS, lat, lon, this.size/2);
        var cst = [], i, j, slen = ortho.length, co = StarMap.CONSTELLATIONS, clen = co.length, halfsize = Math.floor(this.size/2);

        this.constel.remove();
        if (this.eqGrid) { 
            this.eqGrid.remove();
        }
        this.stars.remove();

        for (i = 0; i< clen; ++i) {
            var arc = co[i];
            for (j = arc.length; j--; ) {
                var s = arc[j][0], e = arc[j][1];
                var so = ortho[s], eo = ortho[e];
                if (so[3] || eo[3]) {
                    cst.push('M'+(so[1]+halfsize)+','+(halfsize-so[2])+' L'+(eo[1]+halfsize)+','+(halfsize-eo[2]));
                }
            }
        }

        this.constel.push(this.paper.path(cst.join(' ')).attr({
            'stroke': '#FFF',
            'stroke-width': '1'
        }));

        for (i = 0; i < slen; ++i) {
            var s = ortho[i];
            if (s[3]) {
                this.stars.push(
                    this.paper.circle(s[1]+halfsize, halfsize-s[2], Math.max(3.5-s[0]/2, 0.5)).attr({fill: '#FFF'}));
            }
        }
    }
};
