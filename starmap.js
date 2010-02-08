/**
 * Celestial map component.
 * @constructor
 */
function StarMap (elt, size, prop) {
    this.paper = document.getElementById(elt);
    this.ctx = this.paper.getContext("2d");
    this.prop = prop;

    this.size = size;
    var halfsize = Math.floor(size/2);

    this.drawBg();
}

/**
 * @const
 */
StarMap.DEG2RAD = Math.PI/180.0;

/** @const */
StarMap.STARS=[[2,0.153315,0.321262],[4,0.214373,0.377923],[2,0.277089,0.387169],[3,0.124453,0.349210],[4,0.109105,0.374669],[4,0.152805,0.437325],[3,0.217106,0.451814],[2,0.018302,0.259460],[3,0.086012,0.276017],[4,0.080637,0.303055],[4,0.084312,0.261521],[4,0.103644,0.214995],[4,0.125457,0.207249],[4,-0.047734,0.396615],[3,-0.127389,0.387132],[4,-0.042768,0.407419],[3,-0.048984,0.429200],[5,3.695634,-0.323869],[4,4.869013,-0.277961],[5,3.104510,-0.288432],[5,3.111642,-0.325844],[5,3.624042,-0.300191],[5,2.943970,-0.236295],[5,4.678938,-0.344425],[4,7.197191,-0.335932],[5,6.492389,-0.306284],[5,5.874931,-0.293518],[5,4.425046,-0.258542],[5,-2.446127,-0.790642],[4,-1.567068,-0.819888],[3,-1.472491,-0.822831],[4,-1.407938,-0.802835],[5,-1.917372,-0.746010],[5,-2.011034,-0.745233],[5,-2.518833,-0.742530],[5,-3.565675,-0.792542],[4,-3.214828,-0.854198],[5,-3.115140,-0.840720],[3,-2.607468,-0.824992],[4,-0.768633,0.132259],[2,-0.748742,0.121574],[4,-0.761520,-0.050123],[3,-0.745891,-0.042633],[3,-0.682478,0.027187],[4,-0.655348,0.064482],[2,-0.618034,0.092884],[0,-0.604478,0.077545],[3,-0.591064,0.055967],[3,-0.599455,0.008775],[3,-0.544920,-0.007168],[3,-0.446074,-0.083056],[4,-0.389896,-0.099563],[2,-0.335667,-0.048655],[4,-0.186791,-0.001025],[4,-0.209647,0.012021],[4,-0.228953,-0.068027],[2,-0.254470,-0.002790],[4,-0.252953,-0.121630],[3,-0.217905,-0.012106],[3,-0.148079,-0.066242],[4,-0.154826,-0.119176],[3,-0.143545,-0.138946],[4,-0.096527,-0.079472],[3,-0.110739,-0.186896],[4,-0.109709,-0.198526],[4,-0.116850,-0.210214],[4,-0.099982,-0.052835],[4,-0.074210,-0.182108],[3,-0.080963,-0.177231],[3,-0.971474,-0.467281],[4,-1.112270,-0.429691],[4,-1.089430,-0.460113],[2,-1.131081,-0.464992],[5,-1.071688,-0.485942],[4,-1.305671,-0.500332],[2,-1.164147,-0.526459],[3,-1.311381,-0.531600],[3,-1.163668,-0.535942],[3,-1.365170,-0.566249],[3,-1.134784,-0.585335],[4,0.252876,0.169977],[5,0.298116,0.187242],[2,0.284794,0.207658],[4,0.412200,0.188412],[3,0.388836,0.242484],[4,0.372521,0.246615],[4,0.383634,0.260920],[4,0.444261,0.173868],[2,0.255454,0.183606],[4,0.263110,0.208875],[2,0.756957,0.297790],[3,0.790159,0.376218],[0,0.826867,0.424454],[3,0.997936,0.512659],[2,0.998786,0.336658],[1,0.997944,0.413676],[4,0.999716,0.423825],[2,0.774172,0.402233],[3,0.775953,0.374639],[2,-3.912991,0.161944],[4,-4.106968,0.138741],[-1,-3.279509,0.168980],[2,-2.657170,0.240763],[3,-2.907714,0.271425],[4,-2.726207,0.120378],[4,-3.050561,0.486122],[4,-3.336258,0.485469],[4,-3.261096,0.425384],[4,-4.195054,0.153527],[3,-2.902598,0.347340],[3,-2.385507,0.367834],[3,-2.200616,0.299203],[4,-2.090789,0.338259],[5,0.670816,-0.413742],[4,0.702035,-0.382501],[5,0.706920,-0.335997],[4,0.782715,-0.319942],[4,1.303617,0.795115],[4,0.746905,0.653668],[4,1.085810,0.691396],[4,0.547257,0.643543],[4,0.549636,0.717676],[4,0.779244,0.582507],[4,0.490621,0.576655],[4,0.757964,0.506804],[4,-0.527105,-0.109591],[3,-0.517764,-0.129712],[4,-0.496240,-0.156717],[4,-0.450217,-0.224175],[4,-0.435275,-0.239332],[4,-0.399094,-0.151528],[4,-0.358429,-0.147974],[3,-0.315078,-0.146438],[3,-0.347585,-0.198108],[4,-0.342563,-0.192635],[4,-0.402998,-0.175022],[4,-0.396108,-0.221748],[2,-0.298492,-0.141673],[3,4.235937,-0.700702],[1,2.625613,-0.696536],[4,5.938745,-0.583962],[3,5.145955,-0.597128],[2,5.893326,-0.629665],[3,1.704274,-0.498387],[1,1.946973,-0.571659],[3,2.587637,-0.565392],[2,2.694135,-0.568948],[3,4.378664,-0.592956],[-1,1.110376,-0.495271],[3,1.180002,-0.395888],[2,1.409371,-0.335541],[1,1.818366,-0.438310],[3,4.908485,-0.562775],[4,5.368602,-0.549274],[3,6.843602,-0.564084],[3,8.878472,-0.565484],[2,0.189441,0.580091],[3,0.254885,0.620872],[2,0.020026,0.567498],[2,0.088606,0.537738],[2,0.124353,0.585719],[3,11.726633,-0.514935],[2,-54.833725,-0.474006],[3,18.907471,-0.613037],[3,-39.330620,-0.491720],[3,-16.326486,-0.468759],[2,-11.010328,-0.455301],[2,-4.515992,-0.503673],[0,-3.611286,-0.581699],[2,-3.882782,-0.437807],[-1,-2.754986,-0.587110],[4,-3.676701,-0.420398],[3,-3.775508,-0.412207],[3,-3.789143,-0.384873],[3,-4.101517,-0.388616],[3,-4.105862,-0.380741],[3,-4.968118,-0.358125],[2,-5.628427,-0.331805],[4,-4.259899,-0.296630],[4,-4.017252,-0.296160],[4,-3.966210,-0.286063],[4,-4.108201,-0.310037],[2,-3.525604,-0.328493],[4,-3.158216,-0.343208],[4,-3.099835,-0.359157],[4,-3.544350,-0.375673],[2,-2.833601,-0.385444],[3,-2.426757,-0.384907],[4,-2.711347,-0.342314],[4,-2.680742,-0.316965],[2,-0.367482,0.607836],[4,-0.306907,0.563241],[4,-0.233090,0.543448],[3,-0.242722,0.556605],[3,-0.200791,0.559056],[3,-0.342717,0.707526],[3,-0.154629,0.651898],[4,-0.302278,0.590466],[3,-0.045089,0.804487],[3,-0.452331,0.598948],[4,-0.494280,0.612730],[3,0.042411,-0.077155],[3,0.185391,-0.071534],[3,0.248077,-0.090435],[4,0.362694,0.002867],[4,0.267960,-0.186041],[2,0.095386,-0.158265],[3,0.231025,-0.139984],[3,0.313728,-0.025989],[4,0.334979,0.073961],[3,0.372145,0.028245],[2,0.420053,0.035704],[4,0.376228,0.088492],[4,0.413485,0.077889],[5,-95.990192,-0.788607],[4,1211.464866,-0.812996],[4,5.360806,-0.818603],[4,6.121706,-0.847166],[4,-24.972592,-0.828922],[5,3.030077,-0.853196],[5,2.826885,-0.850870],[5,2.160792,-0.823799],[5,20.140295,-0.779834],[5,3.331676,-0.792250],[4,1.927565,-0.802361],[4,1.905998,-0.794292],[3,-2.701419,-0.636767],[4,-2.175238,-0.563483],[4,-2.103881,-0.569474],[-2,1.219308,-0.146919],[3,1.174390,-0.169638],[1,1.104297,-0.157989],[3,1.269243,-0.214234],[1,1.295162,-0.258358],[4,1.149633,-0.207256],[3,1.245092,-0.291554],[3,1.321188,-0.211036],[3,1.092807,-0.268541],[3,1.313391,-0.248723],[1,1.353840,-0.234486],[4,1.325590,-0.137279],[4,1.280728,-0.149934],[4,1.269576,-0.105444],[3,1.394287,-0.237982],[2,1.455655,-0.261441],[0,1.564420,0.045628],[2,1.476646,0.072465],[4,2.185259,0.189572],[4,2.228958,0.256384],[3,2.203018,0.159763],[3,1.885834,0.080331],[4,2.391846,0.103849],[5,1.843800,0.155235],[4,0.989308,-0.318019],[3,0.996283,-0.392049],[4,1.074959,-0.316647],[3,1.101461,-0.300360],[3,0.881663,-0.319819],[2,0.914920,-0.306438],[3,0.961313,-0.322686],[5,-6.498209,0.154177],[4,-6.325061,0.248197],[4,-16.995781,0.251816],[5,-16.964229,0.238469],[4,-17.342053,0.242556],[5,-18.838648,0.231774],[4,7.567052,-0.161059],[3,6.461938,-0.142264],[3,11.244030,-0.129687],[4,9.446621,-0.201869],[4,13.026495,-0.155557],[4,12.926105,-0.095051],[4,19.639972,-0.085749],[5,115.020419,-0.150796],[4,30.071968,-0.161524],[5,-0.806748,-0.400781],[4,-0.756585,-0.384818],[5,-0.825136,-0.347489],[4,-0.771761,-0.335638],[4,-0.738785,-0.368886],[4,-0.745314,-0.335209],[4,-0.733135,-0.357469],[4,-0.735004,-0.343395],[1,-14.685035,-0.544236],[1,-17.213510,-0.613991],[1,-9.570593,-0.573732],[2,-30.254132,-0.562883],[4,-54.474833,-0.219213],[3,-45.263103,-0.199998],[5,-22.701074,-0.195970],[2,-13.304753,-0.207059],[5,-22.277338,-0.196335],[2,-28.986725,-0.154289],[2,-15.326220,-0.145131],[4,-14.269346,-0.142286],[2,-8.140188,0.347437],[4,-13.560016,0.377444],[3,-0.588144,0.316098],[4,-0.638983,0.269383],[3,-0.665903,0.248954],[2,-0.514403,0.366508],[2,-0.621915,0.415551],[3,-0.709667,0.502604],[3,-0.669110,0.484816],[3,-0.538355,0.432131],[1,-0.462510,0.417082],[2,-0.449911,0.305447],[3,-0.381521,0.270072],[3,-0.421456,0.375548],[3,-0.376892,0.344772],[4,-0.370366,0.357999],[3,-0.401663,0.403293],[5,-0.472115,0.099617],[4,-0.484460,0.098961],[5,-0.469681,0.116721],[5,-0.468620,0.088245],[3,-0.472831,0.128061],[4,-0.478829,0.128760],[3,-0.467268,0.139757],[4,-0.448740,0.141645],[4,-0.457154,0.132313],[3,0.680865,-0.521066],[4,0.624950,-0.482199],[4,0.935674,-0.646131],[5,1.050332,-0.644325],[3,0.891065,-0.606692],[4,0.786604,-0.548309],[3,-13.665024,0.697458],[3,16.007850,0.691543],[3,-3.594055,0.629439],[3,-2.085642,0.565382],[4,-1.715688,0.560775],[2,-1.539926,0.595104],[3,-1.252753,0.645873],[4,-0.913292,0.717749],[2,-1.138109,0.490992],[4,-1.128970,0.522488],[3,-1.028644,0.541519],[2,-1.014917,0.482223],[3,-0.724694,0.670221],[3,-0.912099,0.736358],[3,-0.612280,0.703695],[4,-0.377664,0.087549],[4,-0.388020,0.088646],[3,-0.374324,0.045827],[0,0.216470,-0.545632],[3,0.258520,-0.483514],[3,0.306947,-0.482474],[4,0.332134,-0.442135],[4,0.363476,-0.392818],[4,0.365618,-0.362543],[3,0.409777,-0.366984],[4,0.374191,-0.121532],[3,0.405113,-0.077806],[4,0.376630,-0.163509],[4,0.420339,-0.209133],[4,0.455287,-0.077118],[4,0.528453,-0.337641],[3,0.465024,-0.192187],[3,0.501117,-0.082727],[2,0.794901,-0.044416],[3,0.694306,-0.125477],[4,0.503462,-0.191056],[3,0.529607,-0.085407],[4,0.539706,-0.205722],[4,0.537717,-0.106001],[2,0.571633,-0.118433],[4,0.559196,-0.218149],[4,0.544647,-0.340627],[4,0.577132,-0.212706],[4,0.622662,-0.066882],[4,0.612395,-0.059739],[4,0.679312,-0.265767],[3,0.688303,-0.029264],[3,0.630632,-0.303807],[3,0.649521,-0.305892],[3,0.685836,-0.273214],[4,0.718247,-0.028410],[5,0.390570,-0.321797],[5,0.334664,-0.303929],[5,0.270862,-0.267964],[5,0.486535,-0.321848],[3,0.445416,-0.258496],[5,0.388635,-0.248787],[4,0.386594,-0.290582],[4,0.278476,-0.261384],[3,1.067113,0.198972],[4,1.018140,0.205846],[2,1.105557,0.199035],[1,1.179742,0.144094],[3,1.327697,0.181462],[2,1.212728,0.222895],[3,1.428966,0.194219],[4,1.539289,0.239117],[1,1.610597,0.249570],[3,1.603810,0.216189],[3,1.466815,0.247456],[3,1.415578,0.145351],[3,1.220076,0.113012],[1,1.529616,0.285691],[4,1.490344,0.284710],[4,1.370965,0.270244],[4,1.135092,0.178237],[3,1.261620,0.305359],[3,-0.282197,-0.338140],[4,-0.253700,-0.359463],[2,-0.170330,-0.433616],[1,-0.248787,-0.434408],[3,-0.157143,-0.480375],[3,-0.108723,-0.416737],[4,-0.129699,-0.495905],[4,-0.116411,-0.399167],[3,-0.200568,-0.398914],[4,-1.657964,0.413549],[4,-1.467984,0.388245],[3,-1.571651,0.427705],[2,-1.495052,0.189763],[3,-1.555270,0.168715],[3,-1.301519,0.276632],[2,-1.419649,0.283000],[3,-1.409110,0.353354],[3,-1.093869,0.424540],[3,-1.218250,0.332744],[3,-1.220411,0.126243],[3,-1.218331,0.220222],[3,-1.060903,0.246739],[3,-1.009799,0.260925],[4,-1.136582,0.231884],[4,-0.813324,0.160007],[4,-0.818265,0.181247],[3,-0.901596,0.192295],[4,-0.962496,0.183664],[4,-0.993447,0.190720],[3,-0.967616,0.256407],[4,-1.172488,0.336011],[3,-1.016483,0.337028],[4,0.411141,-0.625735],[5,0.423479,-0.574303],[5,0.327087,-0.580985],[5,0.377515,-0.621288],[3,0.618823,-0.386815],[3,4.109504,-0.108230],[1,2.895725,-0.075704],[3,4.812426,-0.147991],[3,16.957925,-0.285402],[4,64.644279,-0.304853],[3,-5.750358,-0.205012],[3,-3.534759,-0.237150],[4,2.116193,0.049814],[4,2.129407,0.029167],[4,2.184503,0.029667],[3,2.347210,0.051931],[3,2.229999,0.056074],[3,2.645813,0.020197],[3,3.168139,-0.009972],[2,0.056240,-0.799063],[3,0.540806,-0.756827],[4,0.362957,-0.677903],[4,0.319497,-0.682914],[4,0.256140,-0.670038],[2,0.265075,-0.595761],[3,-0.472786,-0.437839],[3,-0.427542,-0.559501],[4,-0.272822,-0.520483],[4,-0.364299,-0.503491],[4,-0.219345,0.430013],[4,-0.197805,0.395154],[4,-0.213560,0.490210],[3,-0.195984,0.469312],[4,-0.211378,0.460756],[4,-0.199978,0.442167],[4,-0.175168,0.406834],[4,-0.235748,0.361156],[4,-0.230937,0.341878],[3,3.986073,0.147335],[2,4.509389,0.174902],[1,4.024679,0.104815],[3,4.361419,0.207246],[3,3.509460,0.230932],[4,7.906107,0.177944],[3,9.983442,0.135466],[2,41.890227,0.127854],[2,9.954611,0.181042],[2,3.318744,0.210497],[3,3.199562,0.086541],[4,2.982628,0.203162],[4,2.836788,0.232543],[3,12.679696,0.092143],[4,11.766032,0.052665],[3,5.193640,0.081395],[2,0.870225,-0.183168],[3,0.786433,-0.197743],[3,0.934406,-0.198443],[2,0.887568,-0.156794],[3,0.813154,-0.142370],[3,0.962832,-0.184248],[3,0.944640,-0.130072],[4,1.021991,-0.144854],[3,0.984435,-0.124271],[5,0.839043,-0.107889],[4,0.837520,-0.115497],[4,0.810862,-0.103950],[4,0.814242,-0.113417],[2,-2.556980,-0.140911],[2,-2.181589,-0.082066],[3,-2.354853,-0.224277],[3,-1.967631,-0.129783],[3,-1.951807,-0.250579],[3,-1.934805,-0.265871],[4,4.908759,0.331756],[4,3.989638,0.317648],[3,6.824691,0.307783],[4,4.803219,0.303786],[4,3.037570,0.328759],[2,-2.711894,-0.438847],[4,-2.544506,-0.399724],[2,-2.436258,-0.395260],[3,-2.122623,-0.327443],[4,-3.028596,-0.418098],[4,-3.029558,-0.416478],[3,-3.186050,-0.425069],[4,-2.787396,-0.460221],[4,-2.891600,-0.471171],[4,-2.339939,-0.435342],[3,-2.246894,-0.452965],[3,-2.242283,-0.488804],[4,-2.288319,-0.417076],[4,-2.162545,-0.443922],[3,-2.112155,-0.411042],[3,-2.127859,-0.370382],[5,-2.087765,-0.361109],[2,-1.971729,-0.375546],[4,-1.941064,-0.389556],[4,-2.106508,-0.333223],[3,-1.730989,-0.348206],[4,-1.923372,-0.309665],[3,-1.813756,-0.302177],[5,-1.759498,-0.305409],[4,-1.893830,-0.312524],[4,1.089509,0.565897],[4,1.287306,0.559138],[4,1.473621,0.457955],[4,1.950340,0.395807],[3,2.726083,0.332679],[3,2.423780,0.381690],[3,2.767278,0.309480],[0,-0.850546,0.351994],[5,-0.823091,0.360712],[4,-0.821506,0.340476],[3,-0.802297,0.299660],[3,-0.770996,0.293267],[4,-0.786557,0.333614],[4,-0.712081,0.345632],[4,-0.916886,0.325541],[4,-0.783624,0.403479],[4,-0.720697,0.355554],[5,0.750775,-0.766453],[5,0.776786,-0.717441],[5,1.045709,-0.763908],[5,0.884265,-0.786090],[5,-0.437378,-0.343481],[6,-0.449653,-0.326080],[5,-0.449588,-0.356075],[5,-0.443957,-0.403910],[5,-0.430486,-0.362093],[5,-0.406650,-0.350505],[5,-0.397885,-0.377730],[5,-0.397918,-0.289974],[6,-0.434864,-0.297896],[4,-0.410915,-0.289181],[5,-0.465431,-0.300318],[4,-0.440092,-0.303629],[5,-0.398908,-0.269118],[5,-0.462607,-0.282955],[4,1.134360,-0.061449],[3,1.067012,-0.054812],[4,1.375529,-0.004300],[4,1.809567,-0.026045],[4,1.729740,-0.032122],[3,1.579147,-0.083542],[4,1.234126,0.021053],[4,1.154957,0.064080],[4,1.109484,0.040101],[4,1.196932,0.086570],[5,12.484482,-0.636518],[5,16.546329,-0.661411],[5,22.333432,-0.641961],[3,31.835179,-0.658485],[4,56.258337,-0.639601],[6,-71.783224,-0.645805],[5,-98.676994,-0.678695],[5,-174.390778,-0.689753],[4,45.571150,-0.703148],[5,20.637220,-0.669697],[4,-26.077663,-0.674011],[5,-9.177490,-0.726360],[2,-12.300138,-0.689024],[3,-14.094408,-0.728318],[5,-5.323003,-0.765772],[3,-7.315605,-0.720544],[5,-5.277516,-0.708399],[3,-9.870171,-0.675881],[4,-6.036558,-0.673169],[4,-5.879839,-0.659173],[4,-1.704331,-0.458149],[4,-1.570895,-0.467961],[4,-1.676766,-0.415986],[4,-1.516686,-0.440583],[4,-3.012362,-0.895168],[3,-0.311764,-0.801007],[4,-0.162737,-0.859857],[2,-1.447160,-0.092478],[4,-1.488589,-0.145998],[4,-1.517851,-0.162469],[4,-1.539105,-0.176664],[4,-1.481608,-0.189552],[2,-1.613266,-0.032251],[3,-1.582442,-0.040972],[4,-1.512256,-0.073186],[3,-1.490169,0.017314],[2,-1.243869,-0.138091],[3,-1.317043,0.081995],[2,-1.074841,0.039877],[4,-1.186469,-0.186357],[4,-1.158665,-0.214155],[3,-1.181195,-0.221689],[5,-1.216625,-0.236402],[4,-1.153643,-0.266705],[2,-1.115815,0.110047],[3,-1.054271,0.023629],[3,-0.997189,0.025589],[4,-0.976476,0.021815],[3,0.732693,0.060824],[4,0.735284,0.077826],[5,0.749395,0.100043],[4,0.754826,0.118485],[3,0.737286,0.048951],[3,0.747590,0.021301],[4,0.762307,0.014960],[2,0.919379,-0.016955],[1,0.901247,-0.010489],[6,0.911466,-0.022642],[6,0.897938,-0.038637],[4,0.897984,-0.042247],[5,0.897544,-0.047068],[2,0.898168,-0.051620],[2,0.884751,-0.002610],[3,0.830240,-0.059800],[0,0.818994,-0.071695],[1,0.858359,0.055468],[0,0.979151,0.064727],[2,0.947956,-0.084585],[4,0.880010,0.051953],[4,1.017252,0.177571],[4,0.975780,0.178812],[4,0.895753,0.083000],[4,1.010453,0.084389],[4,1.033595,0.129596],[4,1.053503,0.124635],[4,-0.827878,-0.718944],[3,-0.575630,-0.738753],[3,-0.552238,-0.651667],[4,-0.777953,-0.664821],[4,-0.794664,-0.603089],[4,-0.903463,-0.594865],[4,-0.963246,-0.620852],[3,-1.064270,-0.633688],[1,-0.505004,-0.539965],[4,-0.348131,-0.641571],[3,-0.453203,-0.651929],[2,0.028885,0.133282],[2,-0.121098,0.133475],[3,-0.251619,0.224855],[4,-0.304192,0.227607],[4,-0.304510,0.152574],[4,-0.358823,0.174567],[4,-0.116135,0.082297],[5,-0.105685,0.076243],[3,-0.161669,0.208597],[2,-0.123283,0.250095],[4,-0.161308,0.106628],[3,-0.173041,0.094804],[2,-0.305285,0.086389],[3,-0.244236,0.054138],[2,-0.169582,0.270020],[4,-0.244729,0.297907],[3,-0.153907,0.218050],[3,0.390629,0.530540],[3,0.399613,0.495996],[2,0.426527,0.504111],[4,0.437582,0.462206],[1,0.477841,0.464832],[2,0.435249,0.373442],[3,0.427507,0.352550],[3,0.438700,0.412754],[4,0.510872,0.447246],[3,0.528704,0.443008],[4,0.494705,0.445179],[3,0.535053,0.389668],[4,0.602827,0.442225],[2,0.571122,0.364071],[2,0.560404,0.285646],[3,0.532600,0.289470],[4,0.574343,0.322905],[4,0.621532,0.449516],[4,0.596661,0.470046],[5,0.567526,0.473715],[4,0.374382,0.458134],[4,0.230087,0.473644],[4,0.390348,0.347440],[3,0.145178,-0.431895],[3,0.201752,-0.456493],[4,0.090405,-0.425350],[3,0.057229,-0.400781],[2,0.057404,-0.386932],[3,0.195206,-0.397123],[3,0.150309,-0.523296],[4,0.253147,-0.427593],[3,0.020535,-0.421876],[3,1.235944,-0.600164],[3,0.946004,-0.477686],[5,1.031274,-0.602699],[5,0.836753,-0.472763],[5,1.137650,-0.541293],[5,0.942585,-0.430639],[5,0.877587,-0.435619],[4,-0.303473,-0.296458],[4,-0.248439,-0.296105],[4,-0.174850,-0.240481],[5,-0.296838,-0.276369],[4,-0.195501,-0.290016],[4,-0.148276,-0.295032],[4,-0.140655,-0.291849],[1,-0.136868,-0.264418],[4,0.138190,0.068962],[4,0.106611,0.066288],[4,0.157624,0.268788],[4,0.162296,0.217885],[4,0.175127,0.242518],[3,0.202278,0.134723],[4,0.224973,0.047924],[4,0.234072,0.080087],[4,0.272740,0.024121],[3,-0.093724,0.028650],[5,-0.086733,0.046996],[4,-0.069996,0.055723],[4,-0.072265,0.010957],[4,-0.043770,0.049139],[4,-0.039187,0.015534],[4,-0.001501,0.059965],[2,1.763766,-0.364003],[3,1.598846,-0.258196],[3,1.642248,-0.220410],[4,1.556912,-0.225036],[4,1.535314,-0.252755],[2,1.799816,-0.215332],[4,1.632544,-0.230291],[3,2.145769,-0.318259],[3,2.189128,-0.297983],[4,2.279842,-0.246641],[3,0.532269,-0.634704],[3,0.620101,-0.606503],[4,0.573705,-0.593759],[4,0.626331,-0.569256],[4,-0.024163,-0.250535],[4,0.128561,-0.261948],[2,-1.729145,-0.200016],[2,-1.685559,-0.174577],[4,-1.631892,-0.171477],[2,-1.560741,-0.227128],[4,-1.629429,-0.248645],[3,-1.759558,-0.260612],[2,-1.742115,-0.231917],[0,-1.500807,-0.234842],[2,-1.455802,-0.251332],[2,-1.362809,-0.308531],[3,-1.352224,-0.344791],[3,-1.335614,-0.387486],[3,-1.234052,-0.396323],[1,-1.104201,-0.393888],[1,-1.122334,-0.335603],[2,-1.136458,-0.337468],[3,-1.045259,-0.335015],[3,-1.055692,-0.365225],[2,-1.079484,-0.354413],[4,-1.107951,-0.350540],[4,-0.880075,-0.127799],[3,-0.857078,-0.072068],[4,-0.812768,-0.041455],[4,-0.830682,-0.079162],[3,-1.975360,0.092231],[2,-1.878073,0.056132],[3,-1.815083,0.039095],[3,-1.859258,0.135400],[3,-1.763422,0.137531],[4,-1.834701,0.159651],[4,-1.905207,0.173361],[3,-1.826348,-0.029943],[3,-1.102910,-0.135192],[4,-1.187389,-0.112582],[4,-1.084567,-0.112833],[3,-1.004255,-0.085498],[3,-0.911086,-0.025302],[3,-0.968433,0.083655],[4,-0.780513,0.036699],[4,4.008484,-0.003243],[5,5.044109,-0.005558],[5,4.793957,-0.061685],[5,3.502021,-0.070847],[5,4.997606,-0.023908],[3,-0.614637,0.163166],[4,-0.636773,0.158508],[4,-0.633859,0.153701],[3,-0.580972,0.171761],[3,-0.925901,-0.332284],[1,-0.899726,-0.309405],[2,-0.912345,-0.266341],[2,-0.974972,-0.271920],[2,-0.884841,-0.225558],[3,-0.818284,-0.239994],[2,-0.783871,-0.233596],[2,-0.758313,-0.266829],[3,-0.743545,-0.246276],[2,-0.734024,-0.185552],[3,-0.694780,-0.157018],[4,-0.691663,-0.408692],[3,-0.687644,-0.370071],[4,-0.694608,-0.140140],[4,-0.591216,-0.382546],[4,-0.578116,-0.317952],[3,-0.941679,-0.185870],[3,-0.775224,-0.186302],[3,0.485493,0.085139],[3,0.479155,0.078955],[3,0.579330,0.109431],[3,0.636436,0.137227],[4,0.694236,0.109613],[4,0.623459,0.077755],[4,0.495516,0.113375],[3,0.586569,0.052312],[4,0.511927,0.003505],[4,0.686173,0.088903],[3,0.646109,0.154294],[4,0.654038,0.157740],[3,0.663691,0.140204],[3,0.663822,0.168960],[0,0.687024,0.145074],[4,0.707531,0.203061],[3,0.906923,0.186625],[5,0.659665,0.203419],[4,0.656591,0.201758],[2,0.541501,0.213514],[1,0.862764,0.254966],[5,-0.712103,-0.418989],[5,-0.745608,-0.491419],[4,-0.651839,-0.446266],[3,-0.888726,-0.424148],[4,-0.952162,-0.424005],[4,-0.881495,-0.456472],[5,-0.793087,-0.488894],[4,-0.772672,-0.497912],[5,-0.612740,-0.535774],[2,-2.157889,-0.683170],[1,-1.372191,-0.687637],[4,-1.955012,-0.653346],[2,-1.775241,-0.617981],[3,0.251836,0.264014],[3,0.290392,0.315176],[4,0.308870,0.304273],[3,-0.093141,-0.557002],[2,-0.225126,-0.580376],[4,0.068929,-0.612297],[4,-1.816275e-4,-0.644174],[4,0.043817,-0.635534],[2,7.839004,0.535998],[1,8.104641,0.597905],[3,-29.701313,0.543322],[3,2.978495,0.613544],[4,10.930906,0.282303],[3,4.621686,0.378860],[3,4.379252,0.393048],[3,2.469268,0.436439],[3,11.008993,0.297107],[3,32.847620,0.442924],[3,9.069375,0.409096],[2,74.284487,0.506174],[1,-4.183790,0.459032],[1,-8.444487,0.531258],[2,-5.400445,0.519738],[3,3.007352,0.484250],[3,2.402444,0.445664],[3,2.030758,0.585735],[4,3.490480,0.510234],[3,3.458632,0.566217],[2,-2.559832,0.755684],[3,-2.135640,0.724330],[4,-1.588703,0.777852],[4,-1.880149,0.806818],[4,-1.389290,0.869857],[4,-1.129226,0.942128],[2,0.343868,0.987239],[2,6.205925,-0.460160],[3,3.632437,-0.515782],[2,2.538594,-0.398276],[1,2.203253,-0.517336],[3,2.960751,-0.368589],[3,4.277652,-0.385085],[3,5.482360,-0.447592],[2,2.787397,-0.520687],[3,-23.011274,-0.005820],[3,49.253479,0.015401],[3,-10.972244,-0.012649],[4,-6.501826,-0.048373],[0,-5.318207,-0.097710],[3,-8.203027,0.029657],[2,-7.326729,0.095929],[4,-3.351852,-0.089895],[4,-3.270471,-0.052412],[3,-4.771472,-0.005199],[4,-3.679132,0.013478],[3,-2.691434,-0.049418],[3,-2.635145,0.016519],[3,-2.052254,0.259603],[2,-1.976573,0.237447],[4,-1.995539,0.280702],[3,-1.893229,0.233586],[4,-1.826584,0.231496],[4,-1.753290,0.238950],[4,-1.719523,0.266557],[4,32.406105,0.057041],[4,-88.002873,0.076358],[3,1.356046,-0.706715],[4,1.803398,-0.682374],[3,1.407331,-0.673965],[4,2.451139,-0.654333],[5,1.919107,-0.720095],[3,1.583525,-0.734655],[3,1.981131,-0.651107],[4,-0.672273,0.218629],[4,-0.596528,0.213283]];

/** @const */
StarMap.CONSTELLATIONS=[[0,1],[2,1],[0,3],[4,3],[4,5],[6,5],[7,8],[9,8],[9,0],[8,0],[8,10],[11,10],[11,12],[9,13],[14,13],[15,13],[15,16],[17,18],[19,20],[20,21],[19,22],[17,23],[23,24],[23,25],[23,26],[21,17],[27,18],[27,21],[28,29],[29,30],[29,31],[30,31],[28,32],[32,33],[33,34],[34,28],[28,35],[35,36],[36,37],[37,38],[38,28],[39,40],[41,42],[40,43],[42,43],[43,44],[44,45],[45,46],[46,47],[48,43],[49,48],[50,51],[51,52],[53,54],[55,56],[57,52],[56,54],[52,56],[56,58],[58,53],[55,59],[60,61],[59,60],[61,62],[63,62],[63,64],[65,64],[66,59],[62,66],[62,59],[67,68],[62,68],[69,70],[71,70],[71,72],[72,73],[73,74],[74,75],[75,76],[76,77],[77,78],[78,79],[79,69],[80,81],[82,81],[81,83],[84,83],[84,85],[84,86],[85,86],[87,83],[88,82],[82,85],[89,88],[89,82],[90,91],[91,92],[92,93],[94,90],[95,92],[95,94],[93,96],[96,95],[92,97],[97,98],[98,91],[99,100],[101,99],[102,103],[104,102],[105,106],[106,107],[108,100],[101,104],[107,105],[103,109],[109,107],[109,110],[102,101],[110,111],[111,102],[111,112],[112,110],[113,114],[114,115],[115,116],[114,116],[117,118],[117,119],[120,121],[120,122],[120,123],[121,118],[118,122],[122,124],[125,126],[127,128],[128,129],[130,126],[127,126],[131,130],[132,131],[133,130],[133,134],[131,134],[135,130],[135,127],[135,136],[137,132],[138,139],[140,141],[141,142],[142,138],[143,144],[144,145],[146,147],[148,139],[149,148],[149,150],[143,151],[152,147],[152,153],[140,154],[154,155],[153,154],[156,157],[158,159],[159,160],[160,156],[161,162],[163,164],[164,165],[162,165],[166,165],[167,166],[167,168],[169,167],[167,170],[169,166],[169,171],[172,171],[172,173],[169,174],[174,175],[175,173],[175,176],[176,177],[177,178],[179,178],[180,178],[178,181],[175,182],[183,182],[183,184],[184,185],[185,173],[173,186],[186,187],[187,188],[188,189],[190,191],[192,191],[192,193],[194,193],[195,190],[195,196],[193,194],[197,193],[194,196],[196,198],[198,195],[190,199],[200,199],[201,202],[202,203],[204,205],[206,205],[206,207],[203,208],[208,204],[209,210],[204,210],[210,211],[212,209],[213,212],[211,213],[201,206],[214,215],[216,215],[216,217],[217,215],[217,218],[217,219],[219,220],[220,221],[221,216],[222,216],[223,221],[221,224],[224,225],[226,227],[226,228],[229,230],[231,230],[232,230],[232,233],[234,230],[233,235],[229,236],[233,237],[238,233],[236,239],[240,241],[240,242],[242,241],[229,241],[239,238],[239,243],[244,243],[245,246],[247,248],[249,250],[249,247],[251,249],[252,250],[252,247],[253,254],[255,256],[255,253],[257,258],[258,259],[259,253],[259,254],[253,256],[260,261],[261,262],[261,263],[261,264],[261,265],[266,267],[266,268],[269,270],[269,266],[268,271],[271,272],[273,272],[270,268],[274,270],[273,274],[275,276],[277,278],[276,279],[280,278],[279,281],[282,280],[281,282],[283,284],[285,286],[287,288],[289,290],[288,291],[292,288],[293,292],[291,293],[294,293],[295,296],[297,298],[298,299],[300,301],[302,301],[302,303],[304,303],[304,305],[300,297],[300,306],[307,306],[307,308],[309,310],[308,310],[305,308],[305,311],[305,300],[312,313],[312,314],[315,313],[315,312],[316,312],[316,317],[316,318],[318,319],[320,314],[319,320],[321,322],[323,324],[324,325],[326,325],[326,321],[326,321],[325,321],[323,325],[327,328],[329,327],[330,329],[331,330],[332,331],[333,332],[333,334],[335,336],[336,337],[337,338],[338,335],[334,339],[334,340],[339,337],[339,341],[342,343],[344,342],[345,346],[346,347],[347,348],[348,349],[349,350],[350,351],[352,353],[354,352],[354,355],[353,356],[351,357],[355,358],[356,359],[360,361],[358,362],[359,363],[362,364],[363,365],[365,366],[364,367],[368,357],[367,369],[366,370],[371,370],[369,372],[371,373],[374,368],[375,374],[372,376],[376,375],[373,377],[377,360],[378,379],[380,379],[378,381],[381,382],[382,383],[383,384],[384,378],[384,385],[385,380],[386,387],[388,386],[389,390],[391,388],[390,392],[393,394],[393,395],[393,396],[393,392],[397,392],[397,398],[390,392],[390,389],[399,400],[401,396],[402,391],[401,403],[401,391],[400,401],[404,405],[406,407],[406,408],[406,409],[406,410],[409,411],[407,412],[407,405],[411,412],[413,414],[415,413],[416,417],[418,416],[419,420],[421,414],[421,415],[420,422],[418,419],[423,424],[425,418],[426,422],[424,427],[428,429],[430,429],[422,418],[430,431],[432,431],[432,433],[422,434],[435,421],[434,435],[427,425],[425,426],[426,433],[436,437],[437,438],[438,439],[439,436],[437,440],[441,442],[443,441],[267,443],[444,269],[445,444],[446,445],[447,446],[448,449],[449,450],[450,451],[452,448],[451,452],[453,451],[442,454],[454,453],[455,456],[457,456],[457,458],[458,459],[459,460],[461,462],[463,462],[464,462],[464,463],[465,466],[467,468],[467,469],[469,470],[470,465],[468,470],[471,466],[471,470],[466,472],[472,473],[474,475],[476,474],[477,478],[475,477],[479,480],[480,474],[481,482],[478,483],[474,484],[479,475],[479,482],[474,483],[483,485],[485,486],[486,478],[480,487],[487,488],[489,480],[490,491],[490,492],[493,494],[493,490],[492,495],[496,493],[495,497],[497,498],[498,496],[499,500],[501,502],[494,500],[494,502],[494,491],[503,504],[505,503],[506,507],[506,503],[504,506],[507,508],[509,510],[509,511],[512,510],[512,511],[510,513],[514,515],[515,516],[516,517],[514,518],[518,519],[514,520],[514,521],[521,522],[521,523],[523,524],[524,525],[523,526],[526,527],[527,528],[526,529],[529,530],[530,531],[531,532],[530,533],[533,517],[534,517],[535,517],[535,536],[537,538],[538,517],[536,537],[539,540],[540,541],[541,542],[543,544],[544,542],[545,543],[546,547],[548,549],[549,550],[550,551],[550,552],[553,546],[547,554],[555,552],[555,554],[549,553],[556,557],[558,557],[559,556],[558,559],[560,561],[561,562],[563,562],[563,564],[564,565],[565,566],[560,564],[565,567],[568,569],[568,570],[570,571],[571,568],[567,569],[567,572],[572,573],[574,575],[576,574],[577,578],[578,579],[578,576],[580,576],[580,581],[580,582],[582,581],[583,581],[584,585],[585,586],[587,586],[587,588],[589,590],[590,591],[591,592],[590,593],[594,595],[596,597],[598,599],[598,600],[600,599],[601,599],[600,602],[602,603],[604,605],[606,604],[605,607],[607,606],[608,609],[609,610],[610,608],[611,612],[613,612],[613,614],[615,614],[616,617],[618,617],[619,616],[611,620],[621,619],[620,622],[620,623],[624,623],[624,625],[626,625],[627,626],[622,621],[611,621],[628,621],[622,628],[629,622],[629,630],[630,631],[632,633],[633,634],[634,635],[636,632],[637,636],[638,637],[639,640],[639,641],[641,642],[642,643],[643,644],[644,645],[646,640],[647,648],[649,632],[649,646],[646,647],[639,650],[651,639],[650,652],[649,652],[653,654],[655,652],[650,656],[656,657],[657,658],[657,654],[658,653],[659,660],[660,661],[661,662],[663,662],[664,665],[663,664],[663,665],[666,665],[662,659],[667,661],[667,668],[669,661],[669,668],[7,670],[670,671],[672,673],[672,674],[673,675],[671,676],[671,677],[678,671],[670,679],[672,673],[680,672],[681,682],[681,683],[684,685],[678,672],[680,681],[686,678],[679,7],[686,684],[671,680],[671,679],[687,688],[687,689],[688,690],[689,691],[692,693],[690,694],[694,692],[691,695],[696,695],[697,698],[697,691],[696,699],[698,700],[701,702],[700,703],[703,701],[699,704],[704,705],[706,705],[691,690],[707,690],[707,708],[688,689],[700,692],[693,709],[710,711],[710,712],[713,714],[712,715],[713,712],[716,712],[710,717],[716,710],[715,710],[711,717],[714,718],[713,718],[719,720],[721,722],[723,724],[720,724],[724,725],[725,722],[722,720],[726,727],[728,727],[726,729],[729,727],[727,730],[730,731],[731,732],[732,733],[733,728],[734,735],[736,737],[737,738],[738,736],[739,737],[740,734],[741,739],[742,740],[742,741],[743,744],[744,745],[746,743],[745,747],[747,748],[748,746],[749,747],[735,749],[151,750],[150,751],[752,753],[753,754],[754,751],[755,752],[756,752],[756,751],[755,750],[750,757],[757,758],[758,759],[760,761],[762,760],[761,763],[763,762],[762,761],[763,760],[764,765],[766,767],[766,767],[767,768],[769,768],[770,769],[770,771],[771,772],[769,773],[773,774],[774,775],[775,776],[776,777],[777,778],[778,779],[780,781],[780,782],[779,783],[784,785],[781,785],[783,784],[786,787],[787,788],[789,786],[789,788],[790,791],[791,792],[793,790],[793,794],[795,794],[795,796],[796,793],[797,616],[792,797],[798,799],[800,798],[801,800],[802,801],[631,802],[631,803],[804,803],[805,806],[807,808],[807,809],[806,809],[808,805],[810,811],[810,812],[813,810],[814,815],[816,817],[815,816],[818,816],[819,816],[820,819],[821,819],[821,822],[822,820],[823,824],[825,826],[826,821],[824,827],[828,829],[826,829],[822,829],[818,830],[830,817],[820,831],[831,823],[823,820],[817,814],[832,833],[834,832],[835,834],[835,836],[835,837],[832,838],[832,839],[833,840],[836,841],[842,835],[843,842],[844,835],[845,843],[846,844],[847,848],[847,849],[850,845],[851,849],[852,847],[848,846],[853,854],[854,855],[854,856],[856,857],[857,858],[858,859],[859,860],[860,861],[861,854],[862,863],[864,862],[865,864],[863,865],[866,867],[866,868],[867,868],[869,870],[869,871],[872,869],[873,872],[871,873],[874,875],[875,876],[875,877],[878,879],[879,880],[879,881],[882,878],[883,884],[884,882],[885,874],[883,886],[876,885],[887,876],[888,887],[886,888],[880,889],[881,890],[877,891],[890,891],[889,892],[893,877],[892,893],[894,895],[895,896],[897,894],[896,897],[898,897],[899,898],[900,899],[901,902],[151,903],[904,151],[903,905],[905,906],[906,907],[907,901],[908,904],[902,908],[909,910],[911,909],[911,912],[913,912],[914,911],[915,914],[913,916],[917,916],[918,914],[919,918],[917,920],[921,919],[917,919],[913,918],[922,923],[924,922],[923,925],[925,926],[926,927],[927,928],[910,929],[929,930],[930,911],[931,932],[933,932],[934,932],[935,932],[936,931],[932,936],[937,932],[934,937],[145,146],[938,939]];

function stereographicProjectPoints(arr, lam1, phi1, rad) {
    function sinSum(cosa, sina, cosb, sinb) {
        return cosa*sinb+sina*cosb;
    }

    function cosSum(cosa, sina, cosb, sinb) {
        return cosa*cosb-sina*sinb;
    }
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

StarMap.prototype.drawBg = function () {
    var size = this.size;
    var halfsize = Math.floor(size/2);
    var ctx = this.ctx;

    ctx.clearRect(0, 0, size, size);
    ctx.fillStyle='#FFF';
    ctx.fillRect(0,0,size, size);
    ctx.beginPath();
    ctx.fillStyle = (this.prop.circleFill || "#000010");
    ctx.arc(halfsize, halfsize, halfsize, 0, 2*Math.PI, true);
    ctx.fill();
}    

StarMap.prototype.setPos = function (lat, lon, time) {
    if (typeof time === 'undefined') {
        time = (new Date()).getTime();
    } else if (typeof time !== 'number') {
        time = time.getTime();
    }
    
    function mjd2jct(mjd) {
        return (mjd - 51544.5) / 36525.0;
    };
    function mod(x, y) {
        return x - y*Math.floor(x/y);
    }
    function gmst(mjd) {
        /** @const */
        var SECS = 86400; // 24*60*60 -- number of seconds in day;
        var mjd0 = Math.floor(mjd), ut = SECS * (mjd - mjd0), t0, t, gmst;
        t0 = mjd2jct(mjd0);
        t = mjd2jct(mjd);
        gmst = 24110.54841 + 8640184.812866 * t0 + 1.0027379093 * ut +
	    (0.093104 - (6.2e-6) * t) * t * t;
        return 2 * Math.PI / SECS * mod(gmst, SECS);
    }

    var mjd = time / 86400000.0 + 40587.0;
    var gms_t = gmst(mjd);

    /** @const */
    var DEG2RAD = StarMap.DEG2RAD;
    lat *= DEG2RAD;
    lon *= DEG2RAD;

    lat += gms_t;

    var ortho = stereographicProjectPoints(StarMap.STARS, lat, lon, this.size/2);
    var cst = [], i, j, slen = ortho.length, co = StarMap.CONSTELLATIONS, clen = co.length, halfsize = Math.floor(this.size/2);
    
    this.drawBg();

    var ctx = this.ctx;
    ctx.beginPath();
    ctx.strokeStyle = 'rgba(255,255,255,0.3)';
    for (j = clen; j--; ) {
        var s = co[j][0], e = co[j][1];
        var so = ortho[s], eo = ortho[e];
        if (so[3] || eo[3]) {
            ctx.moveTo((so[1]+halfsize), (halfsize-so[2]));
            ctx.lineTo((eo[1]+halfsize), (halfsize-eo[2]));
        }
    }
    
    ctx.stroke();
//     this.constel.push(this.paper.path(cst.join(' ')).attr({
//         'stroke': '#FFF',
//         'stroke-opacity': 0.3,
//         'stroke-width': '1'
//     }));
    
    ctx.fillStyle = '#FFF';
    for (i = 0; i < slen; ++i) {
        var s = ortho[i];
        if (s[3]) {
            ctx.beginPath();
            ctx.arc(s[1]+halfsize, halfsize-s[2],
                    Math.max(3.5-s[0]/2, 0.5),
                    0, 2*Math.PI, true);
            ctx.fill();
            //this.paper.circle().attr({fill: '#FFF', 'stroke-width':0}));
        }
    }
};

window['StarMap']=StarMap;
StarMap.prototype['setPos'] = StarMap.prototype.setPos;
