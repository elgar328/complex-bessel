//! Region 2 uniform asymptotic expansion parameter computation.
//!
//! Translation of Fortran ZUNHJ from TOMS 644 / SLATEC (zbsubs.f lines 5017-5731).
//! Computes parameters PHI, ARG, ZETA1, ZETA2, ASUM, BSUM for the uniform
//! asymptotic expansion of J, Y, and H functions in the region |arg(z)| > pi/3.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]
#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::machine::BesselFloat;
use crate::types::SumOption;
use crate::utils::{zabs, zdiv};

// ── Coefficient tables (Fortran DATA blocks, zbsubs.f lines 5068-5442) ──

/// AR(1..14) coefficients (Fortran lines 5068-5076)
#[rustfmt::skip]
const AR: [f64; 14] = [
    1.00000000000000000e+00, 1.04166666666666667e-01,
    8.35503472222222222e-02, 1.28226574556327160e-01,
    2.91849026464140464e-01, 8.81627267443757652e-01,
    3.32140828186276754e+00, 1.49957629868625547e+01,
    7.89230130115865181e+01, 4.74451538868264323e+02,
    3.20749009089066193e+03, 2.40865496408740049e+04,
    1.98923119169509794e+05, 1.79190200777534383e+06,
];

/// BR(1..14) coefficients (Fortran lines 5077-5085)
#[rustfmt::skip]
const BR: [f64; 14] = [
    1.00000000000000000e+00, -1.45833333333333333e-01,
   -9.87413194444444444e-02, -1.43312053915895062e-01,
   -3.17227202678413548e-01, -9.42429147957120249e-01,
   -3.51120304082635426e+00, -1.57272636203680451e+01,
   -8.22814390971859444e+01, -4.92355370523670524e+02,
   -3.31621856854797251e+03, -2.48276742452085896e+04,
   -2.04526587315129788e+05, -1.83844491706820990e+06,
];

/// C(1..105) coefficients (Fortran lines 5086-5152)
#[rustfmt::skip]
const C_COEFFS: [f64; 105] = [
    1.00000000000000000e+00, -2.08333333333333333e-01,
    1.25000000000000000e-01,  3.34201388888888889e-01,
   -4.01041666666666667e-01,  7.03125000000000000e-02,
   -1.02581259645061728e+00,  1.84646267361111111e+00,
   -8.91210937500000000e-01,  7.32421875000000000e-02,
    4.66958442342624743e+00, -1.12070026162229938e+01,
    8.78912353515625000e+00, -2.36408691406250000e+00,
    1.12152099609375000e-01, -2.82120725582002449e+01,
    8.46362176746007346e+01, -9.18182415432400174e+01,
    4.25349987453884549e+01, -7.36879435947963170e+00,
    2.27108001708984375e-01,  2.12570130039217123e+02,
   -7.65252468141181642e+02,  1.05999045252799988e+03,
   -6.99579627376132541e+02,  2.18190511744211590e+02,
   -2.64914304869515555e+01,  5.72501420974731445e-01,
   -1.91945766231840700e+03,  8.06172218173730938e+03,
   -1.35865500064341374e+04,  1.16553933368645332e+04,
   -5.30564697861340311e+03,  1.20090291321635246e+03,
   -1.08090919788394656e+02,  1.72772750258445740e+00,
    2.02042913309661486e+04, -9.69805983886375135e+04,
    1.92547001232531532e+05, -2.03400177280415534e+05,
    1.22200464983017460e+05, -4.11926549688975513e+04,
    7.10951430248936372e+03, -4.93915304773088012e+02,
    6.07404200127348304e+00, -2.42919187900551333e+05,
    1.31176361466297720e+06, -2.99801591853810675e+06,
    3.76327129765640400e+06, -2.81356322658653411e+06,
    1.26836527332162478e+06, -3.31645172484563578e+05,
    4.52187689813627263e+04, -2.49983048181120962e+03,
    2.43805296995560639e+01,  3.28446985307203782e+06,
   -1.97068191184322269e+07,  5.09526024926646422e+07,
   -7.41051482115326577e+07,  6.63445122747290267e+07,
   -3.75671766607633513e+07,  1.32887671664218183e+07,
   -2.78561812808645469e+06,  3.08186404612662398e+05,
   -1.38860897537170405e+04,  1.10017140269246738e+02,
   -4.93292536645099620e+07,  3.25573074185765749e+08,
   -9.39462359681578403e+08,  1.55359689957058006e+09,
   -1.62108055210833708e+09,  1.10684281682301447e+09,
   -4.95889784275030309e+08,  1.42062907797533095e+08,
   -2.44740627257387285e+07,  2.24376817792244943e+06,
   -8.40054336030240853e+04,  5.51335896122020586e+02,
    8.14789096118312115e+08, -5.86648149205184723e+09,
    1.86882075092958249e+10, -3.46320433881587779e+10,
    4.12801855797539740e+10, -3.30265997498007231e+10,
    1.79542137311556001e+10, -6.56329379261928433e+09,
    1.55927986487925751e+09, -2.25105661889415278e+08,
    1.73951075539781645e+07, -5.49842327572288687e+05,
    3.03809051092238427e+03, -1.46792612476956167e+10,
    1.14498237732025810e+11, -3.99096175224466498e+11,
    8.19218669548577329e+11, -1.09837515608122331e+12,
    1.00815810686538209e+12, -6.45364869245376503e+11,
    2.87900649906150589e+11, -8.78670721780232657e+10,
    1.76347306068349694e+10, -2.16716498322379509e+09,
    1.43157876718888981e+08, -3.87183344257261262e+06,
    1.82577554742931747e+04,
];

/// ALFA(1..180) coefficients (Fortran lines 5153-5276)
#[rustfmt::skip]
const ALFA: [f64; 180] = [
   -4.44444444444444444e-03, -9.22077922077922078e-04,
   -8.84892884892884893e-05,  1.65927687832449737e-04,
    2.46691372741792910e-04,  2.65995589346254780e-04,
    2.61824297061500945e-04,  2.48730437344655609e-04,
    2.32721040083232098e-04,  2.16362485712365082e-04,
    2.00738858762752355e-04,  1.86267636637545172e-04,
    1.73060775917876493e-04,  1.61091705929015752e-04,
    1.50274774160908134e-04,  1.40503497391269794e-04,
    1.31668816545922806e-04,  1.23667445598253261e-04,
    1.16405271474737902e-04,  1.09798298372713369e-04,
    1.03772410422992823e-04,  9.82626078369363448e-05,
    9.32120517249503256e-05,  8.85710852478711718e-05,
    8.42963105715700223e-05,  8.03497548407791151e-05,
    7.66981345359207388e-05,  7.33122157481777809e-05,
    7.01662625163141333e-05,  6.72375633790160292e-05,
    6.93735541354588974e-04,  2.32241745182921654e-04,
   -1.41986273556691197e-05, -1.16444931672048640e-04,
   -1.50803558053048762e-04, -1.55121924918096223e-04,
   -1.46809756646465549e-04, -1.33815503867491367e-04,
   -1.19744975684254051e-04, -1.06184319207974020e-04,
   -9.37699549891194492e-05, -8.26923045588193274e-05,
   -7.29374348155221211e-05, -6.44042357721016283e-05,
   -5.69611566009369048e-05, -5.04731044303561628e-05,
   -4.48134868008882786e-05, -3.98688727717598864e-05,
   -3.55400532972042498e-05, -3.17414256609022480e-05,
   -2.83996793904174811e-05, -2.54522720634870566e-05,
   -2.28459297164724555e-05, -2.05352753106480604e-05,
   -1.84816217627666085e-05, -1.66519330021393806e-05,
   -1.50179412980119482e-05, -1.35554031379040526e-05,
   -1.22434746473858131e-05, -1.10641884811308169e-05,
   -3.54211971457743841e-04, -1.56161263945159416e-04,
    3.04465503594936410e-05,  1.30198655773242693e-04,
    1.67471106699712269e-04,  1.70222587683592569e-04,
    1.56501427608594704e-04,  1.36339170977445120e-04,
    1.14886692029825128e-04,  9.45869093034688111e-05,
    7.64498419250898258e-05,  6.07570334965197354e-05,
    4.74394299290508799e-05,  3.62757512005344297e-05,
    2.69939714979224901e-05,  1.93210938247939253e-05,
    1.30056674793963203e-05,  7.82620866744496661e-06,
    3.59257485819351583e-06,  1.44040049814251817e-07,
   -2.65396769697939116e-06, -4.91346867098485910e-06,
   -6.72739296091248287e-06, -8.17269379678657923e-06,
   -9.31304715093561232e-06, -1.02011418798016441e-05,
   -1.08805962510592880e-05, -1.13875481509603555e-05,
   -1.17519675674556414e-05, -1.19987364870944141e-05,
    3.78194199201772914e-04,  2.02471952761816167e-04,
   -6.37938506318862408e-05, -2.38598230603005903e-04,
   -3.10916256027361568e-04, -3.13680115247576316e-04,
   -2.78950273791323387e-04, -2.28564082619141374e-04,
   -1.75245280340846749e-04, -1.25544063060690348e-04,
   -8.22982872820208365e-05, -4.62860730588116458e-05,
   -1.72334302366962267e-05,  5.60690482304602267e-06,
    2.31395443148286800e-05,  3.62642745856793957e-05,
    4.58006124490188752e-05,  5.24595294959114050e-05,
    5.68396208545815266e-05,  5.94349820393104052e-05,
    6.06478527578421742e-05,  6.08023907788436497e-05,
    6.01577894539460388e-05,  5.89199657344698500e-05,
    5.72515823777593053e-05,  5.52804375585852577e-05,
    5.31063773802880170e-05,  5.08069302012325706e-05,
    4.84418647620094842e-05,  4.60568581607475370e-05,
   -6.91141397288294174e-04, -4.29976633058871912e-04,
    1.83067735980039018e-04,  6.60088147542014144e-04,
    8.75964969951185931e-04,  8.77335235958235514e-04,
    7.49369585378990637e-04,  5.63832329756980918e-04,
    3.68059319971443156e-04,  1.88464535514455599e-04,
    3.70663057664904149e-05, -8.28520220232137023e-05,
   -1.72751952869172998e-04, -2.36314873605872983e-04,
   -2.77966150694906658e-04, -3.02079514155456919e-04,
   -3.12594712643820127e-04, -3.12872558758067163e-04,
   -3.05678038466324377e-04, -2.93226470614557331e-04,
   -2.77255655582934777e-04, -2.59103928467031709e-04,
   -2.39784014396480342e-04, -2.20048260045422848e-04,
   -2.00443911094971498e-04, -1.81358692210970687e-04,
   -1.63057674478657464e-04, -1.45712672175205844e-04,
   -1.29425421983924587e-04, -1.14245691942445952e-04,
    1.92821964248775885e-03,  1.35592576302022234e-03,
   -7.17858090421302995e-04, -2.58084802575270346e-03,
   -3.49271130826168475e-03, -3.46986299340960628e-03,
   -2.82285233351310182e-03, -1.88103076404891354e-03,
   -8.89531718383947600e-04,  3.87912102631035228e-06,
    7.28688540119691412e-04,  1.26566373053457758e-03,
    1.62518158372674427e-03,  1.83203153216373172e-03,
    1.91588388990527909e-03,  1.90588846755546138e-03,
    1.82798982421825727e-03,  1.70389506421121530e-03,
    1.55097127171097686e-03,  1.38261421852276159e-03,
    1.20881424230064774e-03,  1.03676532638344962e-03,
    8.71437918068619115e-04,  7.16080155297701002e-04,
    5.72637002558129372e-04,  4.42089819465802277e-04,
    3.24724948503090564e-04,  2.20342042730246599e-04,
    1.28412898401353882e-04,  4.82005924552095464e-05,
];

/// BETA(1..210) coefficients (Fortran lines 5277-5421)
#[rustfmt::skip]
const BETA: [f64; 210] = [
    1.79988721413553309e-02,  5.59964911064388073e-03,
    2.88501402231132779e-03,  1.80096606761053941e-03,
    1.24753110589199202e-03,  9.22878876572938311e-04,
    7.14430421727287357e-04,  5.71787281789704872e-04,
    4.69431007606481533e-04,  3.93232835462916638e-04,
    3.34818889318297664e-04,  2.88952148495751517e-04,
    2.52211615549573284e-04,  2.22280580798883327e-04,
    1.97541838033062524e-04,  1.76836855019718004e-04,
    1.59316899661821081e-04,  1.44347930197333986e-04,
    1.31448068119965379e-04,  1.20245444949302884e-04,
    1.10449144504599392e-04,  1.01828770740567258e-04,
    9.41998224204237509e-05,  8.74130545753834437e-05,
    8.13466262162801467e-05,  7.59002269646219339e-05,
    7.09906300634153481e-05,  6.65482874842468183e-05,
    6.25146958969275078e-05,  5.88403394426251749e-05,
   -1.49282953213429172e-03, -8.78204709546389328e-04,
   -5.02916549572034614e-04, -2.94822138512746025e-04,
   -1.75463996970782828e-04, -1.04008550460816434e-04,
   -5.96141953046457895e-05, -3.12038929076098340e-05,
   -1.26089735980230047e-05, -2.42892608575730389e-07,
    8.05996165414273571e-06,  1.36507009262147391e-05,
    1.73964125472926261e-05,  1.98672978842133780e-05,
    2.14463263790822639e-05,  2.23954659232456514e-05,
    2.28967783814712629e-05,  2.30785389811177817e-05,
    2.30321976080909144e-05,  2.28236073720348722e-05,
    2.25005881105292418e-05,  2.20981015361991429e-05,
    2.16418427448103905e-05,  2.11507649256220843e-05,
    2.06388749782170737e-05,  2.01165241997081666e-05,
    1.95913450141179244e-05,  1.90689367910436740e-05,
    1.85533719641636667e-05,  1.80475722259674218e-05,
    5.52213076721292790e-04,  4.47932581552384646e-04,
    2.79520653992020589e-04,  1.52468156198446602e-04,
    6.93271105657043598e-05,  1.76258683069991397e-05,
   -1.35744996343269136e-05, -3.17972413350427135e-05,
   -4.18861861696693365e-05, -4.69004889379141029e-05,
   -4.87665447413787352e-05, -4.87010031186735069e-05,
   -4.74755620890086638e-05, -4.55813058138628452e-05,
   -4.33309644511266036e-05, -4.09230193157750364e-05,
   -3.84822638603221274e-05, -3.60857167535410501e-05,
   -3.37793306123367417e-05, -3.15888560772109621e-05,
   -2.95269561750807315e-05, -2.75978914828335759e-05,
   -2.58006174666883713e-05, -2.41308356761280200e-05,
   -2.25823509518346033e-05, -2.11479656768912971e-05,
   -1.98200638885294927e-05, -1.85909870801065077e-05,
   -1.74532699844210224e-05, -1.63997823854497997e-05,
   -4.74617796559959808e-04, -4.77864567147321487e-04,
   -3.20390228067037603e-04, -1.61105016119962282e-04,
   -4.25778101285435204e-05,  3.44571294294967503e-05,
    7.97092684075674924e-05,  1.03138236708272200e-04,
    1.12466775262204158e-04,  1.13103642108481389e-04,
    1.08651634848774268e-04,  1.01437951597661973e-04,
    9.29298396593363896e-05,  8.40293133016089978e-05,
    7.52727991349134062e-05,  6.69632521975730872e-05,
    5.92564547323194704e-05,  5.22169308826975567e-05,
    4.58539485165360646e-05,  4.01445513891486808e-05,
    3.50481730031328081e-05,  3.05157995034346659e-05,
    2.64956119950516039e-05,  2.29363633690998152e-05,
    1.97893056664021636e-05,  1.70091984636412623e-05,
    1.45547428261524004e-05,  1.23886640995878413e-05,
    1.04775876076583236e-05,  8.79179954978479373e-06,
    7.36465810572578444e-04,  8.72790805146193976e-04,
    6.22614862573135066e-04,  2.85998154194304147e-04,
    3.84737672879366102e-06, -1.87906003636971558e-04,
   -2.97603646594554535e-04, -3.45998126832656348e-04,
   -3.53382470916037712e-04, -3.35715635775048757e-04,
   -3.04321124789039809e-04, -2.66722723047612821e-04,
   -2.27654214122819527e-04, -1.89922611854562356e-04,
   -1.55058918599093870e-04, -1.23778240761873630e-04,
   -9.62926147717644187e-05, -7.25178327714425337e-05,
   -5.22070028895633801e-05, -3.50347750511900522e-05,
   -2.06489761035551757e-05, -8.70106096849767054e-06,
    1.13698686675100290e-06,  9.16426474122778849e-06,
    1.56477785428872620e-05,  2.08223629482466847e-05,
    2.48923381004595156e-05,  2.80340509574146325e-05,
    3.03987774629861915e-05,  3.21156731406700616e-05,
   -1.80182191963885708e-03, -2.43402962938042533e-03,
   -1.83422663549856802e-03, -7.62204596354009765e-04,
    2.39079475256927218e-04,  9.49266117176881141e-04,
    1.34467449701540359e-03,  1.48457495259449178e-03,
    1.44732339830617591e-03,  1.30268261285657186e-03,
    1.10351597375642682e-03,  8.86047440419791759e-04,
    6.73073208165665473e-04,  4.77603872856582378e-04,
    3.05991926358789362e-04,  1.60315694594721630e-04,
    4.00749555270613286e-05, -5.66607461635251611e-05,
   -1.32506186772982638e-04, -1.90296187989614057e-04,
   -2.32811450376937408e-04, -2.62628811464668841e-04,
   -2.82050469867598672e-04, -2.93081563192861167e-04,
   -2.97435962176316616e-04, -2.96557334239348078e-04,
   -2.91647363312090861e-04, -2.83696203837734166e-04,
   -2.73512317095673346e-04, -2.61750155806768580e-04,
    6.38585891212050914e-03,  9.62374215806377941e-03,
    7.61878061207001043e-03,  2.83219055545628054e-03,
   -2.09841352012720090e-03, -5.73826764216626498e-03,
   -7.70804244495414620e-03, -8.21011692264844401e-03,
   -7.65824520346905413e-03, -6.47209729391045177e-03,
   -4.99132412004966473e-03, -3.45612289713133280e-03,
   -2.01785580014170775e-03, -7.59430686781961401e-04,
    2.84173631523859138e-04,  1.10891667586337403e-03,
    1.72901493872728771e-03,  2.16812590802684701e-03,
    2.45357710494539735e-03,  2.61281821058334862e-03,
    2.67141039656276912e-03,  2.65203073395980430e-03,
    2.57411652877287315e-03,  2.45389126236094427e-03,
    2.30460058071795494e-03,  2.13684837686712662e-03,
    1.95896528478870911e-03,  1.77737008679454412e-03,
    1.59690280765839059e-03,  1.42111975664438546e-03,
];

/// GAMA(1..30) coefficients (Fortran lines 5422-5442)
#[rustfmt::skip]
const GAMA: [f64; 30] = [
    6.29960524947436582e-01, 2.51984209978974633e-01,
    1.54790300415655846e-01, 1.10713062416159013e-01,
    8.57309395527394825e-02, 6.97161316958684292e-02,
    5.86085671893713576e-02, 5.04698873536310685e-02,
    4.42600580689154809e-02, 3.93720661543509966e-02,
    3.54283195924455368e-02, 3.21818857502098231e-02,
    2.94646240791157679e-02, 2.71581677112934479e-02,
    2.51768272973861779e-02, 2.34570755306078891e-02,
    2.19508390134907203e-02, 2.06210828235646240e-02,
    1.94388240897880846e-02, 1.83810633800683158e-02,
    1.74293213231963172e-02, 1.65685837786612353e-02,
    1.57865285987918445e-02, 1.50729501494095594e-02,
    1.44193250839954639e-02, 1.38184805735341786e-02,
    1.32643378994276568e-02, 1.27517121970498651e-02,
    1.22761545318762767e-02, 1.18338262398482403e-02,
];

// Named constants (Fortran lines 5443-5447)
const EX1: f64 = 3.33333333333333333e-01; // 1/3
const EX2: f64 = 6.66666666666666667e-01; // 2/3
use crate::algo::constants::{HPI, PI};
const THPI: f64 = 4.71238898038468986e+00; // 3*pi/2

/// Output of the ZUNHJ parameter computation.
#[derive(Debug, Clone, Copy)]
pub(crate) struct UnhjOutput<T: BesselFloat> {
    pub phi: Complex<T>,
    pub arg: Complex<T>,
    pub zeta1: Complex<T>,
    pub zeta2: Complex<T>,
    pub asum: Complex<T>,
    pub bsum: Complex<T>,
}

/// Compute parameters for the region 2 uniform asymptotic expansion.
///
/// Equivalent to Fortran ZUNHJ in TOMS 644 (zbsubs.f lines 5017-5731).
///
/// # Parameters
/// - `z`: complex argument
/// - `fnu`: order ν > 0
/// - `ipmtr`: `SumOption::Full` = compute all, `SumOption::SkipSum` = skip asum/bsum
/// - `tol`: tolerance for convergence
pub(crate) fn zunhj<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    ipmtr: SumOption,
    tol: T,
) -> UnhjOutput<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let cone = Complex::from(one);

    let rfnu = one / fnu;

    // ── Overflow test: z/fnu too small (Fortran lines 5449-5464) ──
    let test = T::MACH_TINY * T::from_f64(1.0e3);
    let ac = fnu * test;
    if z.re.abs() <= ac && z.im.abs() <= ac {
        let zeta1 = Complex::new(T::from_f64(2.0) * test.ln().abs() + fnu, zero);
        let zeta2 = Complex::new(fnu, zero);
        return UnhjOutput {
            phi: cone,
            arg: cone,
            zeta1,
            zeta2,
            asum: czero,
            bsum: czero,
        };
    }

    // ── Compute in the fourth quadrant (Fortran lines 5466-5477) ──
    let zb = z * rfnu;
    let rfnu2 = rfnu * rfnu;
    let fn13 = fnu.powf(T::from_f64(EX1)); // fnu^(1/3)
    let fn23 = fn13 * fn13; // fnu^(2/3)
    let rfn13 = one / fn13;

    // W2 = 1 - ZB^2 (Fortran lines 5475-5476)
    let w2 = cone - zb * zb;
    let aw2 = zabs(w2);

    if aw2 <= T::from_f64(0.25) {
        // ── Power series for |W2| <= 0.25 (Fortran lines 5482-5579) ──
        return small_w2_branch(w2, aw2, fnu, fn23, rfn13, rfnu, rfnu2, ipmtr, tol);
    }

    // ── |W2| > 0.25: direct computation (Fortran lines 5584-5731) ──
    large_w2_branch(w2, aw2, zb, fnu, fn23, rfn13, rfnu, rfnu2, ipmtr, tol)
}

/// Small |W2| branch: power series (Fortran lines 5482-5579).
fn small_w2_branch<T: BesselFloat>(
    w2: Complex<T>,
    aw2: T,
    fnu: T,
    fn23: T,
    rfn13: T,
    rfnu: T,
    rfnu2: T,
    ipmtr: SumOption,
    tol: T,
) -> UnhjOutput<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    // Power series summation: p[k] = w2^k (Fortran lines 5482-5497)
    let mut p = [czero; 30];
    let mut ap = [zero; 30];

    p[0] = Complex::from(one);
    let mut suma = Complex::new(T::from_f64(GAMA[0]), zero);
    ap[0] = one;

    let mut kmax: usize = 1;
    if aw2 >= tol {
        for k in 1..30 {
            p[k] = p[k - 1] * w2;
            suma = suma + p[k] * T::from_f64(GAMA[k]);
            ap[k] = ap[k - 1] * aw2;
            if ap[k] < tol {
                kmax = k + 1;
                break;
            }
            if k == 29 {
                kmax = 30;
            }
        }
    }

    // ZETA = W2 * SUMA (Fortran lines 5500-5501)
    let zeta = w2 * suma;

    // ARG = ZETA * FN^(2/3) (Fortran lines 5502-5503)
    let arg = zeta * fn23;

    // ZETA2 = sqrt(W2) * FNU (Fortran lines 5504-5507)
    let zeta2 = w2.sqrt() * fnu;

    // ZETA1 = (1 + EX2*(ZETA*SUMA_sqrt)) * ZETA2 (Fortran lines 5504-5511)
    let za = suma.sqrt();
    let ex2_t = T::from_f64(EX2);
    let prod = zeta * za * ex2_t;
    let zeta1 = Complex::new(one + prod.re, prod.im) * zeta2;

    // PHI = sqrt(2*ZA) * RFN13 (Fortran lines 5512-5516)
    let phi = (za + za).sqrt() * rfn13;

    if ipmtr == SumOption::SkipSum {
        return UnhjOutput {
            phi,
            arg,
            zeta1,
            zeta2,
            asum: czero,
            bsum: czero,
        };
    }

    // ── Sum series for ASUM and BSUM (Fortran lines 5521-5579) ──
    let mut bsum = czero;
    for k in 0..kmax {
        bsum = bsum + p[k] * T::from_f64(BETA[k]);
    }

    let mut asum = czero;
    let mut l1: usize = 0;
    let mut l2: usize = 30;
    let btol = tol * (bsum.re.abs() + bsum.im.abs());
    let mut atol_val = tol;
    let mut pp = one;
    let mut ias = false;
    let mut ibs = false;

    if rfnu2 >= tol {
        for _is in 1..7 {
            atol_val = atol_val / rfnu2;
            pp = pp * rfnu2;

            if !ias {
                let mut suma_a = czero;
                for k in 0..kmax {
                    let m = l1 + k;
                    suma_a = suma_a + p[k] * T::from_f64(ALFA[m]);
                    if ap[k] < atol_val {
                        break;
                    }
                }
                asum = asum + suma_a * pp;
                if pp < tol {
                    ias = true;
                }
            }

            if !ibs {
                let mut sumb_b = czero;
                for k in 0..kmax {
                    let m = l2 + k;
                    sumb_b = sumb_b + p[k] * T::from_f64(BETA[m]);
                    if ap[k] < atol_val {
                        break;
                    }
                }
                bsum = bsum + sumb_b * pp;
                if pp < btol {
                    ibs = true;
                }
            }

            if ias && ibs {
                break;
            }
            l1 += 30;
            l2 += 30;
        }
    }

    // Finalize (Fortran lines 5575-5578)
    asum = asum + one;
    let pp_final = rfnu * rfn13;
    bsum = bsum * pp_final;

    UnhjOutput {
        phi,
        arg,
        zeta1,
        zeta2,
        asum,
        bsum,
    }
}

/// Large |W2| branch: direct computation (Fortran lines 5584-5731).
fn large_w2_branch<T: BesselFloat>(
    w2: Complex<T>,
    aw2: T,
    zb: Complex<T>,
    fnu: T,
    fn23: T,
    rfn13: T,
    rfnu: T,
    rfnu2: T,
    ipmtr: SumOption,
    tol: T,
) -> UnhjOutput<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    // W = sqrt(W2) (Fortran line 5585)
    let w = w2.sqrt();
    let mut wr = w.re;
    let mut wi = w.im;
    // Clamping (Fortran lines 5586-5587)
    if wr < zero {
        wr = zero;
    }
    if wi < zero {
        wi = zero;
    }
    let w_c = Complex::new(wr, wi); // clamped w as Complex

    // ZA = (1+W)/ZB, ZC = log(ZA) (Fortran lines 5588-5594)
    let za = zdiv(w_c + one, zb);
    let zc = za.ln();
    let mut zcr = zc.re;
    let mut zci = zc.im;
    // Clamping (Fortran lines 5592-5594)
    if zci < zero {
        zci = zero;
    }
    if zci > T::from_f64(HPI) {
        zci = T::from_f64(HPI);
    }
    if zcr < zero {
        zcr = zero;
    }

    // ZTH = 1.5*(ZC - W) (Fortran lines 5595-5596)
    let zthr = (zcr - wr) * T::from_f64(1.5);
    let zthi = (zci - wi) * T::from_f64(1.5);

    // ZETA1 = ZC * FNU, ZETA2 = W * FNU (Fortran lines 5597-5600)
    let zeta1 = Complex::new(zcr * fnu, zci * fnu);
    let zeta2 = w_c * fnu;

    // Compute AZTH, ANG for ZTH^(2/3) (Fortran lines 5601-5612)
    let zth = Complex::new(zthr, zthi);
    let azth = zabs(zth);
    let ang = {
        let thpi_t = T::from_f64(THPI);
        let hpi_t = T::from_f64(HPI);
        let gpi_t = T::from_f64(PI);
        if zthr >= zero && zthi < zero {
            thpi_t
        } else if zthr == zero {
            hpi_t
        } else {
            let a = (zthi / zthr).atan();
            if zthr < zero { a + gpi_t } else { a }
        }
    };

    let ex2_t = T::from_f64(EX2);
    let pp = azth.powf(ex2_t);
    let ang_scaled = ang * ex2_t;
    let zetar = pp * ang_scaled.cos();
    let mut zetai = pp * ang_scaled.sin();
    if zetai < zero {
        zetai = zero;
    }

    // ARG = ZETA * FN^(2/3) (Fortran lines 5614-5615)
    let zeta_c = Complex::new(zetar, zetai);
    let arg = zeta_c * fn23;

    // RTZTR, RTZTI = ZTH/ZETA (Fortran line 5616)
    let rzth = zdiv(zth, zeta_c);

    // ZA = RZTH / W (Fortran line 5617)
    let za = zdiv(rzth, w_c);

    // PHI = sqrt(2*ZA) * RFN13 (Fortran lines 5620-5622)
    let phi = (za + za).sqrt() * rfn13;

    if ipmtr == SumOption::SkipSum {
        return UnhjOutput {
            phi,
            arg,
            zeta1,
            zeta2,
            asum: czero,
            bsum: czero,
        };
    }

    // ── Compute ASUM and BSUM for |W2| > 0.25 (Fortran lines 5624-5731) ──
    let raw = one / aw2.sqrt();
    let tfn = w_c.conj() * (raw * raw * rfnu);
    let razth = one / azth;
    let rzth_base = zth.conj() * (razth * razth * rfnu);

    let zc_init = rzth_base * T::from_f64(AR[1]); // AR(2), 0-based [1]

    let raw2 = one / aw2;
    let t2 = w2.conj() * (raw2 * raw2);

    // UP(2) (Fortran lines 5641-5644, 0-based index 1)
    let c2_val = T::from_f64(C_COEFFS[1]); // C(2)
    let c3_val = T::from_f64(C_COEFFS[2]); // C(3)
    let up1 = Complex::new(t2.re * c2_val + c3_val, t2.im * c2_val) * tfn;

    let mut bsum = up1 + zc_init;
    let mut asum = czero;

    if rfnu < tol {
        // Skip refinement
        asum = Complex::from(one);
        let bsum_div = zdiv(-bsum * rfn13, rzth);
        return UnhjOutput {
            phi,
            arg,
            zeta1,
            zeta2,
            asum,
            bsum: bsum_div,
        };
    }

    // Full refinement (Fortran lines 5650-5724)
    let mut przth = rzth_base;
    let mut ptfn = tfn;

    let mut up = [czero; 14];
    up[0] = Complex::from(one); // UP(1) = 1
    up[1] = up1; // UP(2) already computed

    let mut cr = [czero; 14];
    let mut dr = [czero; 14];

    let mut pp = one;
    let btol = tol * (bsum.re.abs() + bsum.im.abs());
    let mut ks: usize = 0;
    let mut kp1: usize = 2; // Fortran KP1 starts at 2
    let mut l: usize = 3; // Fortran L starts at 3 (1-based), maps to 0-based index 2
    let mut ias = false;
    let mut ibs = false;

    // DO 210 LR=2,12,2 (Fortran lines 5663-5724)
    let mut lr = 2;
    while lr <= 12 {
        let lrp1 = lr + 1;

        // Compute two more CR, DR, UP (Fortran DO 160 K=LR,LRP1)
        for _k in lr..=lrp1 {
            ks += 1;
            kp1 += 1;
            l += 1;
            // Horner evaluation of C_COEFFS polynomial in t2
            let mut za = Complex::new(T::from_f64(C_COEFFS[l - 1]), zero);
            for _j in 1..kp1 {
                l += 1;
                za = za * t2 + Complex::new(T::from_f64(C_COEFFS[l - 1]), zero);
            }
            // PTFN = PTFN * TFN (Fortran lines 5681-5683)
            ptfn = ptfn * tfn;
            // UP(KP1) (Fortran lines 5684-5685)
            up[kp1 - 1] = ptfn * za;
            // CR(KS) (Fortran lines 5686-5687)
            cr[ks - 1] = przth * T::from_f64(BR[ks]); // BR(KS+1), 0-based [ks]
            // PRZTH = PRZTH * RZTH (Fortran lines 5688-5690)
            przth = przth * rzth_base;
            // DR(KS) (Fortran lines 5691-5692)
            dr[ks - 1] = przth * T::from_f64(AR[ks + 1]); // AR(KS+2), 0-based [ks+1]
        }

        pp = pp * rfnu2;

        if !ias {
            // SUMA (Fortran lines 5696-5707)
            let mut suma_a = up[lrp1 - 1]; // UP(LRP1)
            let mut ju = lrp1; // Fortran 1-based index
            for cr_k in cr[..lr].iter() {
                ju -= 1;
                suma_a = suma_a + *cr_k * up[ju - 1];
            }
            asum = asum + suma_a;
            let test_a = suma_a.re.abs() + suma_a.im.abs();
            if pp < tol && test_a < tol {
                ias = true;
            }
        }

        if !ibs {
            // SUMB (Fortran lines 5710-5721)
            let mut sumb_b = up[lr + 2 - 1]; // UP(LR+2)
            sumb_b = sumb_b + up[lrp1 - 1] * zc_init;
            let mut ju = lrp1;
            for dr_k in dr[..lr].iter() {
                ju -= 1;
                sumb_b = sumb_b + *dr_k * up[ju - 1];
            }
            bsum = bsum + sumb_b;
            let test_b = sumb_b.re.abs() + sumb_b.im.abs();
            if pp < btol && test_b < btol {
                ibs = true;
            }
        }

        if ias && ibs {
            break;
        }
        lr += 2;
    }

    // Label 220: finalize (Fortran lines 5726-5730)
    asum = asum + one;
    let bsum_div = zdiv(-bsum * rfn13, rzth);

    UnhjOutput {
        phi,
        arg,
        zeta1,
        zeta2,
        asum,
        bsum: bsum_div,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::SumOption;
    use num_complex::Complex64;

    const TOL: f64 = 2.220446049250313e-16;

    #[test]
    fn zunhj_small_w2_real() {
        // z = fnu * 0.99 + 0i → w2 = 1 - 0.99^2 ≈ 0.0199 (small branch)
        let fnu = 100.0;
        let z = Complex64::new(fnu * 0.99, 0.0);
        let result = zunhj(z, fnu, SumOption::Full, TOL);
        assert!(result.phi.re.is_finite());
        assert!(result.asum.re.is_finite());
        assert!(result.bsum.re.is_finite());
        // ASUM should be close to 1 (leading term)
        assert!((result.asum.re - 1.0).abs() < 0.1);
    }

    #[test]
    fn zunhj_large_w2_real() {
        // z = 10 + 0i, fnu = 100 → ZB = 0.1, W2 = 1-0.01 = 0.99 (large branch)
        let z = Complex64::new(10.0, 0.0);
        let result = zunhj(z, 100.0, SumOption::Full, TOL);
        assert!(result.phi.re.is_finite());
        assert!(result.zeta1.re > 0.0);
        assert!(result.zeta2.re > 0.0);
        assert!((result.asum.re - 1.0).abs() < 0.5);
    }

    #[test]
    fn zunhj_complex_argument() {
        let z = Complex64::new(5.0, 3.0);
        let result = zunhj(z, 50.0, SumOption::Full, TOL);
        assert!(result.phi.re.is_finite());
        assert!(result.phi.im.is_finite());
        assert!(result.arg.re.is_finite());
    }

    #[test]
    fn zunhj_ipmtr_1() {
        // ipmtr=1: compute only phi, arg, zeta1, zeta2 (no asum/bsum)
        let z = Complex64::new(5.0, 3.0);
        let result = zunhj(z, 50.0, SumOption::SkipSum, TOL);
        assert!(result.phi.re.is_finite());
        // asum/bsum should be zero
        assert_eq!(result.asum.re, 0.0);
        assert_eq!(result.bsum.re, 0.0);
    }

    #[test]
    fn zunhj_overflow_precheck() {
        // Very small z → overflow pre-check path
        let z = Complex64::new(1e-310, 1e-310);
        let result = zunhj(z, 1.0, SumOption::Full, TOL);
        assert_eq!(result.phi.re, 1.0);
        assert!(result.zeta1.re > result.zeta2.re);
    }

    #[test]
    fn zunhj_imaginary_argument() {
        // Pure imaginary z (used in uni2/unk2 after rotation)
        let z = Complex64::new(0.0, 50.0);
        let result = zunhj(z, 100.0, SumOption::Full, TOL);
        assert!(result.phi.re.is_finite());
        assert!(result.zeta1.re.is_finite());
    }
}
