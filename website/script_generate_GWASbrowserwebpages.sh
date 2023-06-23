# TL - 281221
# LAST MAJ: 210623 # PLEASE ALSO CHANGE INFO SHOWN IN L.62 REGARDING LAST UPDATE

#### 1/ generate the index

# print the frame of all webpages > tmpframe
echo -e "<!DOCTYPE html>
<html lang="en">

<head>
  <meta charset="utf-8">
  <link rel="stylesheet" href="style_menu.css">
  <title>rose GWAS browser</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="author" content="Thibault Leroy">
<link rel="preconnect" href="https\:\/\/fonts.googleapis.com">
<link rel="preconnect" href="https\:\/\/fonts.gstatic.com" crossorigin>
<link href="https\:\/\/fonts.googleapis.com\/css2?family=Amatic+SC:wght@700\&family=Damion\&family=Pacifico\&display=swap" rel="stylesheet"> 
</head>

<body>
  
<div class="container">
  
</div>

<nav>
  <ul>
    <li><a href="index.html">roseGWASbrowser</a></li>
    <li class="menu-deroulant">
      <a href="#">Traits</a>
      <ul class="sous-menu">
        <li><a href="./blackspot/index.html">Black Spot Disease</a></li>
        <li><a href="./petalcolor/index.html">Colors of Petals</a></li>
        <li><a href="./numberofpetals/index.html">Number of Petals & Flowers</a></li>
        <li><a href="./numberofpricklesacicules/index.html">Number of Prickles or Acicules</a></li>
        <li><a href="./plantgrowtharchitecture/index.html">Plant Growth & Architecture</a></li>
        <li><a href="./recurrentblooming/index.html">Recurrent Blooming</a></li>
        <li><a href="./scent/index.html">Scent</a></li>
      </ul>
    </li>
    <li><a href="help.html">Help</a></li>
    <li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>
    <li><a href="mailto:thibault.leroy@inrae.fr">Contact</a></li>
  </ul>
</nav>
<br> <br> 
<script>
</script>
<img src="./IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">
</body>
</html>" > tmpframe

# the tmpframe is the general index
cp tmpframe index.html
rose_to_keep=$(echo "<img src="./IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
project_partners_to_add=$(echo "<img src="./IMG/project_partners_logo.png" alt="projectpartners" class="center50">")
head -n +$(($(grep -n "$rose_to_keep" index.html | cut -d: -f1 | bc)-1)) index.html > index.html2
echo "<h1 style=\"text-align:center; font-size: 200%;\"; > Welcome to the rose GWAS browser</h1wq!>" >> index.html2
echo -e "$rose_to_keep" >> index.html2
echo -e "$project_partners_to_add" >> index.html2
tail -n +$(( 1 + $(grep -n "$rose_to_keep" index.html | cut -d: -f1) )) index.html >> index.html2
echo "<h3 style=\"text-align:center; color:#808080\";>© Thibault Leroy - IRHS - INRAE Univ.Angers InstitutAgro - All rights reserved<br>Last update: June 21, 2023</h3>" >> index.html2
mv index.html2 index.html

### Generate a help webpage
cp tmpframe help.html
old_IMG=$(echo "<img src="./IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
new_PDF=$(echo "<iframe src=\"./PDF/GWASbrowser_english_short.pdf\" style=\"width:100%; height:400pc\">")
head -n +$(($(grep -n "$old_IMG" help.html | cut -d: -f1 | bc)-1)) help.html > help.html2
echo -e "$new_PDF" >> help.html2
tail -n +$(( 1 + $(grep -n "$old_IMG" help.html | cut -d: -f1) )) help.html >> help.html2
mv help.html2 help.html

#### 2/ generate the index for each trait
for i in blackspot petalcolor numberofpetals numberofpricklesacicules plantgrowtharchitecture recurrentblooming scent; do
	mkdir $i
	### generate the index for each trait
	## first here the idea is to use the same frame but to change the different paths with sed
	cp tmpframe $i/index.html
	sed "s/\/$i//g" $i/index.html | sed 's/\/blackspot/.\/blackspot/g' | sed 's/\/numberofpetals/.\/numberofpetals/g' | sed 's/\/numberofpricklesacicules/.\/numberofpricklesacicules/g' | sed 's/\/petalcolor/.\/petalcolor/g' | sed 's/\/recurrentblooming/.\/recurrentblooming/g' | sed 's/\/scent/.\/scent/g' | sed 's/\/plantgrowtharchitecture/.\/plantgrowtharchitecture/g' | sed 's/style_menu.css/..\/style_menu.css/g' | sed 's/.\/IMG/..\/IMG/g' |  sed 's/.\/PDF/..\/PDF/g' | sed 's/index.html/..\/index.html/g' | sed 's/\/\.\.\//\//g' | sed 's/\.\/\.\.\//\.\//g'> $i/index.html2
	mv $i/index.html2 $i/index.html
	## then, the idea is to add the second menu specific for each trait, to do that I first print a menu in a tmpfile for each trait. Then from the index, I copy all before the help menu and after, plus the tmpfile for each trait itself. This becomes the new tmpfile
	# blackspot
	if [ $i == "blackspot" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
	      <a href="#">Phenotypes</a>
	      <ul class="sous-menu">
	        <li><a href="./blackspot3years.html">Black Spot (average visual rating score, 3 years)</a></li>
	        <li><a href="./blackspot2014.html">Black Spot (visual rating score, year: 2014)</a></li>
	        <li><a href="./blackspot2015.html">Black Spot (visual rating score, year: 2015)</a></li>
	        <li><a href="./blackspot2016.html">Black Spot (visual rating score, year: 2016)</a></li>
	      </ul>
    	</li>" > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html
	elif [ $i == "numberofpetals" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
	      <a href="#">Phenotype</a>
  	    <ul class="sous-menu">
  	      <li><a href="./meanpetalsraw.html">Number of Petals - Quantitative</a></li>
        <li><a href="./meanpetalsbinary10.html">Number of Petals - Qualitative (Binary: <10 vs. >10)</a></li>
        <li><a href="./meanpetalsbinary25.html">Number of Petals - Qualitative (Binary <25 vs. >25)</a></li>
        <li><a href="./nbflowersperinflo.html">Number of Flowers per infloresence - Quantitative</a></li>
      </ul>
    </li>"  > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html
	elif [ $i == "numberofpricklesacicules" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
	      <a href="#">Phenotype</a>
  	    <ul class="sous-menu">
  	      <li><a href="./pricklesquantitative.html">Prickles - Quantitative</a></li>
        <li><a href="./pricklesqualitative.html">Prickles - Qualitative (Binary: Few vs. Many)</a></li>
        <li><a href="./aciculesquantitative.html">Acicules - Quantitative</a></li>
        <li><a href="./aciculesqualitative.html">Acicules - Qualitative (Binary: Few vs. Many)</a></li>
      </ul>
    </li>"  > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html
	elif [ $i == "petalcolor" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
	      <a href="#">Phenotype</a>
  	    <ul class="sous-menu">
  	      <li><a href="./whitecol.html">White petals (Binary: White or Not) </a></li>
        <li><a href="./pinkcol.html">Pink petals (Binary: Pink or Not)</a></li>
        <li><a href="./purplecol.html">Purple petals (Binary: Purple or Not)</a></li>
        <li><a href="./pinkpurplered.html"> PPR petals (Binary: Pink-Purple-Red vs. other)</a></li>
        <li><a href="./unicolor.html">Unicolor petals (Binary: Unicolor or Not)</a></li>
      </ul>
    </li>"  > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html

	elif [ $i == "plantgrowtharchitecture" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
      <a href="#">Parameters</a>
      <ul class="sous-menu">
        <li><a href="./meanheight.html">Height (mean, 3 years)</a></li>
        <li><a href="./sdheight.html">Height (SD, 3 years)</a></li>
        <li><a href="./height2012.html">Height (year 2012)</a></li>
        <li><a href="./height2013.html">Height (year 2013</a></li>
        <li><a href="./height2014.html">Height (year 2014)</a></li>
        <li><a href="./meancircumference.html">Circumference (mean, 3 years)</a></li>
        <li><a href="./sdcircumference.html">Circumference (sd, 3 years)</a></li>
        <li><a href="./circumference2012.html">Circumference (year 2012)</a></li>
        <li><a href="./circumference2013.html">Circumference (year 2013)</a></li>
        <li><a href="./circumference2014.html">Circumference (year 2014)</a></li>
        <li><a href="./shrubbyornot.html">Shrubby Architecture (Binary: Shrubby Or Not)</a></li>
        <li><a href="./bushyornot.html">Bushy Architecture (Binary: Bushy Or Not)</a></li>
      </ul>
    </li>"  > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html


	elif [ $i == "recurrentblooming" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
      <a href="#">Parameters</a>
      <ul class="sous-menu">
        <li><a href="./reb_mag.html">nGMM - Reb_Mag</a></li>
        <li><a href="./first_reb.html">nGMM - First_Reb</a></li>
        <li><a href="./cont_flo.html">nGMM - Cont_Flo</a></li>
        <li><a href="./area_peak.html">nGMM - Area_Peak</a></li>
        <li><a href="./max_peak.html">nGMM - Max_Peak</a></li>
        <li><a href="./max_reb.html">nGMM - Max_Reb</a></li>
        <li><a href="./rat_peak.html">nGMM - Rat_Peak</a></li>
        <li><a href="./first_flo.html">nGMM - Fist_Flo</a></li>
        <li><a href="./nb_clust.html">GMM - Nb_Clust</a></li>
        <li><a href="./lmax_p.html">GMM - LMax_P</a></li>
        <li><a href="./rat_p.html">GMM - Rat_P</a></li>
        <li><a href="./difmax_m.html">GMM - DifMax_M</a></li>
        <li><a href="./min_m.html">GMM - Min_M</a></li>
        <li><a href="./min_v.html">GMM - Min_V</a></li>
      </ul>
    </li>"  > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html

	elif [ $i == "scent" ]; then
		echo "Currently working on $i"
		echo -e "      <li class="menu-deroulant">
      <a href="#">Compounds</a>
      <ul class="sous-menu">
        <li><a href="./benzylalcohol2yearquantitative.html">BenzylAlcohol (Quantitative) </a></li>
        <li><a href="./benzylalcohol2yearqualitative0p1.html">BenzylAlcohol (Qualitative, threshold:0.1) </a></li>
        <li><a href="./benzylalcohol2yearqualitative5.html">BenzylAlcohol (Qualitative, threshold:5) </a></li>
        <li><a href="./citronellol2yearquantitative.html">betaCitronellol (Quantitative)</a></li>
        <li><a href="./citronellol2yearqualitative0p1.html">betaCitronellol (Qualitative, threshold:0.1)</a></li>
        <li><a href="./citronellol2yearqualitative15.html">betaCitronellol (Qualitative, threshold:15)</a></li>
        <li><a href="./dihydrobionol2yearquantitative.html">Dihydrobetaionol (Quantitative)</a></li>
        <li><a href="./dihydrobionol2yearqualitative0p1.html">Dihydrobetaionol (Qualitative, threshold:0.1)</a></li>
        <li><a href="./dihydrobionol2yearqualitative10.html">Dihydrobetaionol (Qualitative, threshold:10)</a></li>
        <li><a href="./dmt2yearquantitative.html">DMT (Quantitative) </a></li>
        <li><a href="./dmt2yearqualitative0p1.html">DMT (Qualitative, threshold:0.1) </a></li>
        <li><a href="./e2hexenal2yearquantitative.html">E2hexenal (Quantitative) </a></li>
        <li><a href="./e2hexenal2yearqualitative1.html">E2hexenal (Qualitative, threshold:1) </a></li>
        <li><a href="./e2hexenal2yearqualitative0p1.html">E2hexenal (Qualitative, threshold:0.1) </a></li>
        <li><a href="./ebfarnesene2yearquantitative.html">Ebetafarnesene (Quantitative) </a></li>
        <li><a href="./ebfarnesene2yearqualitative0p1.html">Ebetafarnesene (Qualitative, threshold:0.1) </a></li>
        <li><a href="./elemol2yearquantitative.html">Elemol (Quantitative) </a></li>
        <li><a href="./elemol2yearqualitative0p1.html">Elemol (Qualitative, threshold:0.1) </a></li>
        <li><a href="./elemol2yearqualitative5.html">Elemol (Qualitative, threshold:5) </a></li>
        <li><a href="./eugenol2yearquantitative.html">Eugenol (Quantitative) </a></li>
        <li><a href="./eugenol2yearqualitative0p1.html">Eugenol (Qualitative, threshold:0.1) </a></li>
        <li><a href="./farnesol2yearquantitative.html">eeFarnesol (Quantitative) </a></li>
        <li><a href="./farnesol2yearqualitative0p1.html">eeFarnesol (Qualitative, threshold:0.1) </a></li>
        <li><a href="./geranial2yearquantitative.html">Geranial (Quantitative)</a></li>
        <li><a href="./geranial2yearqualitative0p1.html">Geranial (Qualitative, threshold:0.1)</a></li>
        <li><a href="./geranial2yearqualitative10.html">Geranial (Qualitative, threshold:10)</a></li>
        <li><a href="./geraniol2yearquantitative.html">Geraniol (Quantitative)</a></li>
        <li><a href="./geraniol2yearqualitative10.html">Geraniol (Qualitative, threshold:10)</a></li>
        <li><a href="./geraniol2yearqualitative50.html">Geraniol (Qualitative, threshold:50)</a></li>
        <li><a href="./geranylacetate2yearquantitative.html">GeranylAcetate (Quantitative) </a></li>
        <li><a href="./geranylacetate2yearqualitative0p1.html">GeranylAcetate (Qualitative, threshold:0.1) </a></li>
        <li><a href="./geranylacetate2yearqualitative5.html">GeranylAcetate (Qualitative, threshold:5) </a></li>
        <li><a href="./germacrened2yearquantitative.html">GermacreneD (Quantitative) </a></li>
        <li><a href="./germacrened2yearqualitative0p1.html">GermacreneD (Qualitative, threshold:0.1) </a></li>
        <li><a href="./nerol2yearquantitative.html">Nerol (Quantitative) </a></li>
        <li><a href="./nerol2yearqualitative10.html">Nerol (Qualitative, threshold:10) </a></li>
        <li><a href="./pe2yearquantitative.html">2PhenylEthanol (Quantitative)</a></li>
        <li><a href="./pe2yearqualitative10.html">2PhenylEthanol (Qualitative, threshold:10)</a></li>
        <li><a href="./pe2yearqualitative100.html">2PhenylEthanol (Qualitative, threshold:100)</a></li>
        <li><a href="./tmb2yearquantitative.html">TMB (Quantitative) </a></li>
        <li><a href="./tmb2yearqualitative0p1.html">TMB (Qualitative, threshold:0.1) </a></li>
        <li><a href="./sumfarnesol2yearquantitative.html">SUM FARNESOL & derivates (Quantitative) </a></li>
        <li><a href="./sumgeraniol2yearquantitative.html">SUM GERANIOL & derivates (Quantitative) </a></li>
      </ul>
    </li>"  > tmpsecondmenu
		#echo "$second_menu"
		menu_to_replace=$(echo "<li><a href="help.html">Help</a></li>")
		head -n +$(($(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1 | bc)-1)) $i/index.html > $i/index.html2
		cat tmpsecondmenu >> $i/index.html2
		tail -n +$(( 1 + $(grep -n "$menu_to_replace" $i/index.html | cut -d: -f1) )) $i/index.html >> $i/index.html2
		mv $i/index.html2 $i/index.html
	else
		echo "ERROR: the trait $i is not included in the different traits possible, revise my code!"
	fi
done

### 3 - Generate all webpages per trait
for i in blackspot petalcolor numberofpetals numberofpricklesacicules plantgrowtharchitecture recurrentblooming scent; do
	echo "working on $i"
	# cp index (previous step) of the trait for all webpages
	if [ $i == "blackspot" ]; then
		#pwd
		cd $i
		#pwd
		for jshort in 2014 2015 2016 3years; do
			j=$(echo "blackspot""$jshort")
			jshort2=$(echo "$jshort" | sed 's/years/yrs/g')
			cp index.html $j.html
			# replace IMG by PDF
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/circlize_TN""$jshort2""_generalslidwin1000000.pdf\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" index.html | cut -d: -f1 | bc)-1)) index.html > $j.html2
			echo -e "$new_PDF" >> $j.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" index.html | cut -d: -f1) )) index.html >> $j.html2
			mv $j.html2 $j.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$j.general.html">General</a></li>
        <li><a href="./$j.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$j.additive.html">Additive</a></li>
        <li><a href="./$j.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$j.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$j.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$j.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$j.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>" > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" index.html | cut -d: -f1 | bc)-1)) $j.html > $j.html2
			cat tmpthirdmenu >> $j.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" index.html | cut -d: -f1) )) $j.html >> $j.html2
			mv $j.html2 $j.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				oldpdffile=$(echo "circlize_TN""$jshort2""_generalslidwin1000000.pdf")
				newpdffile=$(echo "circlize_TN""$jshort2""_""$kvaredited""slidwin1000000.pdf")
				sed "s/$oldpdffile/$newpdffile/g" $j.html > $j"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..
	elif [ $i == "numberofpetals" ]; then
		#pwd
		cd $i
		for jshort in meanpetalsraw meanpetalsbinary10 meanpetalsbinary25 nbflowersperinflo; do
			j=$(echo "$jshort" | sed 's/meanpetalsraw/mean_petals/g' | sed 's/meanpetalsbinary10/mean_petals_binary_less10/g' | sed 's/meanpetalsbinary25/mean_petals_binary_more25/g' | sed 's/nbflowersperinflo/nb_flowers_per_inflorescence/g')
			cp index.html $jshort.html
			# replace IMG by PDF
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/circlize_""$j""_generalslidwin1000000.pdf\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			echo -e "$new_PDF" >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$jshort.general.html">General</a></li>
        <li><a href="./$jshort.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$jshort.additive.html">Additive</a></li>
        <li><a href="./$jshort.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$jshort.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$jshort.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>"  > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			cat tmpthirdmenu >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				oldpdffile=$(echo "circlize_""$j""_generalslidwin1000000.pdf")
				newpdffile=$(echo "circlize_""$j""_""$kvaredited""slidwin1000000.pdf")
				sed "s/$oldpdffile/$newpdffile/g" $jshort.html > $jshort"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..
	elif [ $i == "numberofpricklesacicules" ]; then
		#pwd
		cd $i
		for jshort in pricklesquantitative pricklesqualitative aciculesquantitative aciculesqualitative; do
			j=$(echo "$jshort" | sed 's/pricklesquantitative/Nb_prickles_branch/g' | sed 's/pricklesqualitative/Nb_prickles_branch_binaryTatiana/g' | sed 's/aciculesquantitative/Nb_acicules_branch/g' | sed 's/aciculesqualitative/Nb_acicules_branch_binaryTatiana/g')
			cp index.html $jshort.html
			# replace IMG by PDF
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/circlize_""$j""_generalslidwin1000000.pdf\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			echo -e "$new_PDF" >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$jshort.general.html">General</a></li>
        <li><a href="./$jshort.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$jshort.additive.html">Additive</a></li>
        <li><a href="./$jshort.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$jshort.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$jshort.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>"  > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			cat tmpthirdmenu >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				oldpdffile=$(echo "circlize_""$j""_generalslidwin1000000.pdf")
				newpdffile=$(echo "circlize_""$j""_""$kvaredited""slidwin1000000.pdf")
				sed "s/$oldpdffile/$newpdffile/g" $jshort.html > $jshort"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..
	elif [ $i == "petalcolor" ]; then
		#pwd
		cd $i
		for jshort in pinkpurplered whitecol pinkcol purplecol  unicolor; do
			j=$(echo "$jshort" | sed 's/whitecol/main_petal_color_white_vs_other/g' | sed 's/pinkcol/main_petal_color_rose_vs_other/g' | sed 's/purplecol/main_petal_color_purple_vs_other/g' | sed 's/pinkpurplered/main_petal_color_redpinkpurple_vs_other/g' | sed 's/unicolor/color_distribution_UnicolorOrNot/g')
			cp index.html $jshort.html
			# replace IMG by PDF
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/circlize_""$j""_generalslidwin1000000.pdf\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			echo -e "$new_PDF" >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$jshort.general.html">General</a></li>
        <li><a href="./$jshort.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$jshort.additive.html">Additive</a></li>
        <li><a href="./$jshort.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$jshort.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$jshort.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>"  > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			cat tmpthirdmenu >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				oldpdffile=$(echo "circlize_""$j""_generalslidwin1000000.pdf")
				newpdffile=$(echo "circlize_""$j""_""$kvaredited""slidwin1000000.pdf")
				sed "s/$oldpdffile/$newpdffile/g" $jshort.html > $jshort"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..
	elif [ $i == "plantgrowtharchitecture" ]; then
		#pwd
		cd $i
		for jshort in meanheight sdheight height2012 height2013 height2014 meancircumference sdcircumference circumference2012 circumference2013 circumference2014 shrubbyornot bushyornot; do
			j=$(echo "$jshort" | sed 's/meanheight/Height_mean/g' | sed 's/sdheight/Height_sd/g' | sed 's/height2012/Height_2012/g' | sed 's/height2013/Height_2013/g' | sed 's/height2014/Height_2014/g' | sed 's/meancircumference/Circumference_mean/g' | sed 's/sdcircumference/Circumference_sd/g' | sed 's/circumference2012/Circumference_2012/g' | sed 's/circumference2013/Circumference_2013/g' | sed 's/circumference2014/Circumference_2014/g' | sed 's/shrubbyornot/Plant_Type_ArbustifOuAutre/g' | sed 's/bushyornot/Plant_Type_BuissonnantOuAutre/g' )
			cp index.html $jshort.html
			# replace IMG by PDF
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/circlize_""$j""_generalslidwin1000000.pdf\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			echo -e "$new_PDF" >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$jshort.general.html">General</a></li>
        <li><a href="./$jshort.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$jshort.additive.html">Additive</a></li>
        <li><a href="./$jshort.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$jshort.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$jshort.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$jshort.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>"  > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1 | bc)-1)) $jshort.html > $jshort.html2
			cat tmpthirdmenu >> $jshort.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" $jshort.html | cut -d: -f1) )) $jshort.html >> $jshort.html2
			mv $jshort.html2 $jshort.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				oldpdffile=$(echo "circlize_""$j""_generalslidwin1000000.pdf")
				newpdffile=$(echo "circlize_""$j""_""$kvaredited""slidwin1000000.pdf")
				sed "s/$oldpdffile/$newpdffile/g" $jshort.html > $jshort"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..
	elif [ $i == "recurrentblooming" ]; then
		cd $i
		for j in reb_mag first_reb cont_flo area_peak max_peak max_reb rat_peak first_flo nb_clust lmax_p rat_p difmax_m min_m min_v; do
			jvaredited=$(echo "$j" | sed 's/reb_mag/Reb_Mag/g' | sed 's/first_reb/First_Reb/g' | sed 's/cont_flo/Cont_Flo/g' | sed 's/area_peak/Area_Peak_Inflo/g' | sed 's/max_peak/Max_Peak/g' | sed 's/max_reb/Max_Reb/g' | sed 's/rat_peak/Rat_Peak/g' | sed 's/first_flo/First_Flo/g' | sed 's/nb_clust/Nb_Clust/g' | sed 's/lmax_p/Lmax_P/g' | sed 's/rat_p/Rat_P/g' | sed 's/difmax_m/DifMax_M/g' | sed 's/min_m/Min_M/g' | sed 's/min_v/Min_V/g')
			cp index.html $j.html
			# replace IMG by PDF
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/circlize_PF_""$jvaredited""_generalslidwin_1000000.pdf\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" $j.html | cut -d: -f1 | bc)-1)) $j.html > $j.html2
			echo -e "$new_PDF" >> $j.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" $j.html | cut -d: -f1) )) $j.html >> $j.html2
			mv $j.html2 $j.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$j.general.html">General</a></li>
        <li><a href="./$j.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$j.additive.html">Additive</a></li>
        <li><a href="./$j.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$j.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$j.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$j.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$j.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>"  > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" $j.html | cut -d: -f1 | bc)-1)) $j.html > $j.html2
			cat tmpthirdmenu >> $j.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" $j.html | cut -d: -f1) )) $j.html >> $j.html2
			mv $j.html2 $j.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				oldpdffile=$(echo "circlize_PF_""$jvaredited""_generalslidwin_1000000.pdf")
				newpdffile=$(echo "circlize_PF_""$jvaredited""_""$kvaredited""slidwin_1000000.pdf")
				sed "s/$oldpdffile/$newpdffile/g" $j.html > $j"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..	
	elif [ $i == "scent" ]; then
		cd $i
		for j in pe2yearquantitative pe2yearqualitative10 pe2yearqualitative100 geraniol2yearquantitative geraniol2yearqualitative10 geraniol2yearqualitative50 dmt2yearquantitative dmt2yearqualitative0p1 farnesol2yearquantitative farnesol2yearqualitative0p1 citronellol2yearquantitative citronellol2yearqualitative0p1 citronellol2yearqualitative15 geranial2yearquantitative geranial2yearqualitative0p1 geranial2yearqualitative10 geranylacetate2yearquantitative geranylacetate2yearqualitative0p1 geranylacetate2yearqualitative5 germacrened2yearquantitative germacrened2yearqualitative0p1 benzylalcohol2yearquantitative benzylalcohol2yearqualitative0p1 benzylalcohol2yearqualitative5 dihydrobionol2yearquantitative dihydrobionol2yearqualitative0p1 dihydrobionol2yearqualitative10 e2hexenal2yearquantitative e2hexenal2yearqualitative0p1 e2hexenal2yearqualitative1 ebfarnesene2yearquantitative ebfarnesene2yearqualitative0p1 elemol2yearquantitative elemol2yearqualitative0p1 elemol2yearqualitative5 eugenol2yearquantitative eugenol2yearqualitative0p1 nerol2yearquantitative nerol2yearqualitative10 tmb2yearquantitative tmb2yearqualitative0p1 sumfarnesol2yearquantitative sumgeraniol2yearquantitative; do
			cp index.html $j.html
			# replace IMG by PDF
			jvaredited=$(echo "$j" | sed 's/sumfarnesol2yearquantitative/Sum_farnesolfarnesalfarnesylacetate_2years/g' | sed 's/sumgeraniol2yearquantitative/Sum_geraniolgeranialgeranylacetatecitronellolnerol_2years/g' | sed 's/pe2yearquantitative/X2phenylethanol_2years/g' | sed 's/pe2yearqualitative100/X2phenylethanol_binary100_2years/g' | sed 's/pe2yearqualitative10/X2phenylethanol_binary10_2years/g' | sed 's/geraniol2yearquantitative/geraniol_2years/g' | sed 's/geraniol2yearqualitative10/geraniol_binary10_2years/g' | sed 's/geraniol2yearqualitative50/geraniol_binary50_2years/g' | sed 's/dmt2yearquantitative/X35dimethoxytoluene_2years/g' | sed 's/dmt2yearqualitative0p1/X35dimethoxytoluene_binary0p1_2years/g' | sed 's/farnesol2yearquantitative/EEfarnesol_2years/g' | sed 's/farnesol2yearqualitative0p1/EEfarnesol_binary0p1_2years/g' | sed 's/citronellol2yearquantitative/bcitronellol_2years/g' | sed 's/citronellol2yearqualitative0p1/bcitronellol_binary0p1_2years/g'  | sed 's/citronellol2yearqualitative15/bcitronellol_binary15_2years/g'  | sed 's/geranial2yearquantitative/geranial_2years/g' | sed 's/geranial2yearqualitative0p1/geranial_binary0p1_2years/g' | sed 's/geranial2yearqualitative10/geranial_binary10_2years/g' | sed 's/geranylacetate2yearquantitative/geranylacetate_2years/g' | sed 's/geranylacetate2yearqualitative0p1/geranylacetate_binary0p1_2years/g'  | sed 's/geranylacetate2yearqualitative5/geranylacetate_binary5_2years/g'  | sed 's/germacrened2yearquantitative/germacreneD_2years/g' | sed 's/germacrened2yearqualitative0p1/germacreneD_binary0p1_2years/g' | sed 's/benzylalcohol2yearqualitative10/benzylalcohol_binary10_2years/g' | sed 's/benzylalcohol2yearquantitative/benzylalcohol_2years/g' | sed 's/benzylalcohol2yearqualitativeOp1/benzylalcohol_binary0p1_2years/g'  | sed 's/benzylalcohol2yearqualitative5/benzylalcohol_binary5_2years/g' | sed 's/dihydrobionol2yearquantitative/dihydrobionol_2years/g' | sed 's/dihydrobionol2yearqualitative0p1/dihydrobionol_binary0p1_2years/g' | sed 's/dihydrobionol2yearqualitative10/dihydrobionol_binary10_2years/g' | sed 's/e2hexenal2yearquantitative/E2hexenal_2years/g' | sed 's/e2hexenal2yearqualitative0p1/E2hexenal_binary0p1_2years/g' | sed 's/e2hexenal2yearqualitative1/E2hexenal_binary1_2years/g' | sed 's/ebfarnesene2yearquantitative/Ebfarnesene_2years/g' | sed 's/ebfarnesene2yearqualitative0p1/Ebfarnesene_binary0p1_2years/g' | sed 's/elemol2yearquantitative/elemol_2years/g' | sed 's/elemol2yearqualitative0p1/elemol_binary0p1_2years/g' | sed 's/elemol2yearqualitative5/elemol_binary5_2years/g' | sed 's/eugenol2yearquantitative/eugenol_2years/g' | sed 's/eugenol2yearqualitative0p1/eugenol_binary0p1_2years/g' | sed 's/nerol2yearquantitative/nerol_2years/g' | sed 's/nerol2yearqualitative10/nerol_binary10_2years/g' |  sed 's/tmb2yearquantitative/TMB_2years/g' | sed 's/tmb2yearqualitative0p1/TMB_binaryPresAbs_2years/g' )
			newpdffile=$(echo "circlize_""$jvaredited""_generalslidwin1000000.pdf")
			old_IMG=$(echo "<img src="../IMG/francois_foucard_P1140924.png" alt="foucardrosedet" class="center35">")
			new_PDF=$(echo "<iframe src=\"../PDF/""$newpdffile""\" style=\"width:100%; height:400pc\">")
			head -n +$(($(grep -n "$old_IMG" $j.html | cut -d: -f1 | bc)-1)) $j.html > $j.html2
			echo -e "$new_PDF" >> $j.html2
			tail -n +$(( 1 + $(grep -n "$old_IMG" $j.html | cut -d: -f1) )) $j.html >> $j.html2
			mv $j.html2 $j.html
			# replace contact menu by the menu allowing the choice of models 
 			echo -e "    <li class="menu-deroulant">
      <a href="#">Models</a>
      <ul class="sous-menu">
        <li><a href="./$j.general.html">General</a></li>
        <li><a href="./$j.diplo-general.html">Diploidized general</a></li>
        <li><a href="./$j.additive.html">Additive</a></li>
        <li><a href="./$j.diplo-additive.html">Diplodized additive</a></li>
        <li><a href="./$j.simplex-refdom.html">Simplex Dominant (Ref)</a></li>
        <li><a href="./$j.simplex-altdom.html">Simplex Dominant (Alt)</a></li>
        <li><a href="./$j.duplex-refdom.html">Duplex Dominant (Ref)</a></li>
        <li><a href="./$j.duplex-altdom.html">Duplex Dominant (Alt)</a></li>
      </ul>
    </li>" > tmpthirdmenu
			menu_to_replace2=$(echo "<li><a href="https://github.com/roseGWASbrowser/PopGen_GWAS_19thcentury_roses" target="_blank">Github</a></li>")
			head -n +$(($(grep -n "$menu_to_replace2" $j.html | cut -d: -f1 | bc)-1)) $j.html > $j.html2
			cat tmpthirdmenu >> $j.html2
			tail -n +$(( 1 + $(grep -n "$menu_to_replace2" $j.html | cut -d: -f1) )) $j.html >> $j.html2
			mv $j.html2 $j.html
			# then for each 
			for k in general diplo-general additive diplo-additive simplex-refdom simplex-altdom duplex-refdom duplex-altdom; do
				kvaredited=$(echo "$k" | sed 's/simplex-refdom/1\.dom\.ref/g' | sed 's/simplex-altdom/1\.dom\.alt/g' | sed 's/duplex-refdom/2\.dom\.ref/g' | sed 's/duplex-altdom/2\.dom\.alt/g' | sed 's/diplo-general/diplo\.general/g' | sed 's/diplo-additive/diplo\.additive/g')
				newpdffile2=$(echo "circlize_""$jvaredited""_""$kvaredited""slidwin1000000.pdf")
				sed "s/$newpdffile/$newpdffile2/g" $j.html > $j"."$k.html
			done
		done
		rm tmpthirdmenu
		cd ..	
	else
		echo "ERROR: the trait $i is not included in the different traits possible, revise my code!"
	fi
done
