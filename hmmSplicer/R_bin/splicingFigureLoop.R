
#Store some metric for each mutation position that helps to sort them according to how well they support and expression effect for each splice mutation
#Each record corresponds to one mutation that may occur in several libraries.  If the mutation is real and affects expression of the involved junction...
#The average or cumulative difference from the median for each mutated library versus all other libraries could be calculated
#Initialize values
id_data[,"cum_median_diff"]=NA
id_data[,"avg_median_diff"]=NA
id_data[,"max_recip_diff"]=NA
id_data[,"max_gain_canon"]=NA
id_data[,"max_loss_canon"]=NA
id_data[,"max_gain_alter"]=NA
id_data[,"max_loss_alter"]=NA

#Loop to process each mutation and generate SVG files
for(i in 1:length(id_data[,1])){
#for(i in 1:10){
#for(i in 508:508){

  #Get the mutation name and position
  mut_name = id_data[i,"mid"]
  mut_pos = id_data[i,"pos"]

  #Get the gene name(s), chr, strand, donor_acceptor status, and coordinate range of interest
  gene_name = id_data[i,"gene_string"]
  donor_acceptor = id_data[i,"donor_acceptor"]
  chr = id_data[i,"chr"]
  strand = id_data[i,"strand"]

  #Get the coordinate range of interest (the outer coords of all junctions associated with the mutation position
  grand_start = id_data[i,"min_start"]
  grand_end = id_data[i,"max_end"]

  #Get the list of mutatated libraries (by id and name)
  mut_lib_list=strsplit(id_data[i,"mutation_lib_list"],split=",")
  mut_lib_list=mut_lib_list[[1]][1:(length(mut_lib_list[[1]]))]
  mut_lib_name_list=id_data[i,"mutation_lib_name_list"]
  mut_lib_name_array=strsplit(id_data[i,"mutation_lib_name_list"],split=",")
  mut_lib_name_array=mut_lib_name_array[[1]][1:(length(mut_lib_name_array[[1]]))]

  #Get the expression values for each known junction
  #Order the lists according to intron size
  #Assign short names to the junctions
  if (!is.na(id_data[i,"known_junction_ids_string"])){
    kj_list=strsplit(as.character(id_data[i,"known_junction_ids_string"]),split=",")
    kj_list=kj_list[[1]][1:(length(kj_list[[1]]))]

    #Make sure the known junction is actually present in the expression data object
    kj_list=kj_list[which(kj_list %in% exp_data[, "JID"])]

    if (length(which(exp_data[, "JID"] %in% kj_list)) > 0){
      intron_sizes = exp_data[kj_list,"Intron_Size"]
      o = order(intron_sizes)
      kj_list=kj_list[o]
      kj_names = paste("CJ_", seq(1:length(kj_list)), sep="")
      names(kj_list) = kj_names
    }else{
      kj_list=NULL
      kj_names=NULL
    }
  }else{
    kj_list=NULL
    kj_names=NULL
  }

  #Get the expression values for each alternate junction
  if (!is.na(id_data[i,"alternate_junction_ids_string"])){
    aj_list=strsplit(as.character(id_data[i,"alternate_junction_ids_string"]),split=",")
    aj_list=aj_list[[1]][1:(length(aj_list[[1]]))]
    if (length(which(exp_data[, "JID"] %in% aj_list)) > 0){
      intron_sizes = exp_data[aj_list,"Intron_Size"]
      o = order(intron_sizes)
      aj_list=aj_list[o]
      aj_names = paste("AJ_", seq(1:length(aj_list)), sep="")
      names(aj_list) = aj_names
    }else{
      aj_list=NULL
      aj_names=NULL
    }
  }else{
    aj_list=NULL
    aj_names=NULL
  }

  #Combine the lists
  j_list = c(kj_list, aj_list)
  j_names = c(kj_names, aj_names)

  if (length(j_list) > 0){
    #Get the normalized expression values for the list of junctions defined
    tmp=exp_data[j_list,lib_ids]
    j_exp_list=t(tmp[j_list,lib_ids])
    colnames(j_exp_list)=j_names

    #Set the boxplot colors
    bp_colors=c(rep(kj_color, length(kj_list)), rep(aj_color, length(aj_list)))
    
    #Define the plot name
    plot_name = paste(donor_acceptor, " mutation at chr", mut_name," (",strand,")", " (", gene_name, ")", sep="")

    #Set up the plotting space
   
    if (active == 1){
      #Using devSVGTips
      outfile = paste(outdir, chr, "_", mut_pos, ".active.svg", sep="")
      devSVGTips(outfile, toolTipMode=2, toolTipFontSize=12, title=plot_name, sub.special=TRUE, width=12, height=15)
    }else{
      #Using Cairo...
      outfile = paste(outdir, chr, "_", mut_pos, ".svg", sep="")
      Cairo(width=800, height=450*2, file=outfile, type="svg", pointsize=12, units="px", dpi=72)
    }
    par(font.main=2, font.axis=2, font.lab=2, mfrow=c(2,1), mar=(c(10, 4, 4, 2)+0.1))
    
    #First, draw the splicing model diagram
    drawModel(ec_data, jc_data, mut_pos, chr, strand, grand_start, grand_end, j_list, j_names, plot_name, kj_color, aj_color, intron_color, exon_color, mutant_color, mut_lib_name_list, active)

    #Next, create a box plot showing the expression distribution of each junction associated with the mutation site
    #Name the plot with the mutation position and gene name(s)
    #Color the boxplots according to whether they correspond to canonical or alternative junctions
    #Mark the values for the mutated libraries on each boxplot
    xlabel="Junction identity"
    if (length(j_list) == 1){
      xlabel=j_names[1]
    }
    ylabel=paste("Normalized junction expression level (", length(lib_ids), " libraries)", sep="")
    max_y = max(j_exp_list)
    max_yp = max_y*0.1

    bp=boxplot(j_exp_list, ylab=ylabel, xlab=xlabel, main="Expression distribution of each junction depicted above", col=bp_colors, ylim=c(0,max_y+max_yp), font.main=2, font.axis=2, font.lab=2)
    legend_text = paste("Mutant libraries (n = ", length(mut_lib_list), ")", sep="")
    legend("topleft", legend_text, pch=8, col=mutant_color)

    #Draw points for the expression level of each mutant library for each junction
    #points(x=rep(1:length(j_list), length(mut_lib_list)), y=t(exp_data[j_list,mut_lib_list]), col=mutant_color, pch=8, cex=2)
    pch_list=c(8,7,9,10,11,12,13,14,15,16,17,18)
    pch_count = 0
    if (length(mut_lib_list) > 12){pch_list=rep(8,length(mut_lib_list))}
    for (lib in mut_lib_list){
      pch_count = pch_count+1
      pch_select = pch_list[pch_count]
      points(x=1:length(j_list), y=exp_data[j_list,lib], col=mutant_color, pch=pch_select, cex=2)
    }

    #If there are both canonical and alternative junctions, draw a dividing line between them
    if (length(kj_list) > 0 & length(aj_list) > 0){
      div = (length(kj_list)) + 0.5
      abline(v=div, lty=2, col="black")
    }

    #If active was specified, add pop-up name for the library name
    if (active == 1){
      pch_count = 0
      for (l in 1:length(mut_lib_list)){
        lib = mut_lib_list[l]
        lib_name = mut_lib_name_array[l]
        pch_count = pch_count+1
        pch_select = pch_list[pch_count]
        url_text = paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=", chr, ":", mut_pos-1, "-", mut_pos+1, sep="")
        title_text = "Mutant Library"
        desc1_text = paste("Library Name: ", lib_name, sep="")

        for (j in 1:length(j_list)){
          desc2_text = paste("Junction: ", j_list[j], sep="")
          setSVGShapeToolTip(title=title_text, desc1=desc1_text, desc2=desc2_text)
          setSVGShapeURL(url_text, target="_blank") 
          points(x=j, y=exp_data[j_list[j],lib], col="red", pch=1, cex=2)
        }
      }
    }

    dev.off()

    #For each junction in 'j_list', calculate the:
    #Cumulative absolute difference from the median for each mutant library in 'mut_lib_list'
    #Average absolute difference from the median for each mutant library
    #Max gain of a canonical junction
    #Max loss of a canonical junction
    #Max gain of an alternative junction
    #Max loss of an alternative junction
    #If there is loss of a canonical and gain of an alternative junction within a single library, the absolute sum of these
    max_cmd = 0
    max_events = 1
    max_gain_canon = 0
    max_loss_canon = 0
    max_gain_alter = 0
    max_loss_alter = 0
    max_recip_diff = 0
    for (lib in mut_lib_list){
      cmd_lib = 0
      events_lib = 0
      max_gain_canon_lib = 0
      max_loss_canon_lib = 0
      max_gain_alter_lib = 0
      max_loss_alter_lib = 0
      for (jid in j_list){
        mut_lib_exp = exp_data[jid,lib]
        exclude_list=lib_ids[which(!lib_ids %in% lib)]
        if (length(exclude_list) > 0){
          med_j_exp = as.numeric(median(exp_data[jid,exclude_list]))
          md = mut_lib_exp - med_j_exp
          cmd_lib = cmd_lib+abs(md)
          events_lib = events_lib+1
        
          #Max gain and loss of a canonical junction
          if (length(kj_list) > 0){
            if (jid %in% kj_list){
              if (md > max_gain_canon){max_gain_canon = md; max_gain_canon_lib = md}
              if (md < max_loss_canon){max_loss_canon = md; max_loss_canon_lib = md}
            }
          }
          #Max gain and loss of an alternative junction
          if (length(aj_list) > 0){
            if (jid %in% aj_list){
              if (md > max_gain_alter){max_gain_alter = md; max_gain_alter_lib = md}
              if (md < max_loss_alter){max_loss_alter = md; max_loss_alter_lib = md}
            }
          }
        }
        #Summarize within this lib
        if (cmd_lib > max_cmd){max_cmd = cmd_lib}
        if (events_lib > max_events){max_events = events_lib}
        if (max_loss_canon_lib < 0 & max_gain_alter_lib > 0){
          recip = abs(max_loss_canon_lib) + max_gain_alter_lib
          if (recip > max_recip_diff){max_recip_diff = recip}
        }
      }
    }
    id_data[i,"cum_median_diff"] = cmd_lib
    id_data[i,"avg_median_diff"] = cmd_lib/max_events
    id_data[i,"max_gain_canon"] = max_gain_canon
    id_data[i,"max_loss_canon"] = max_loss_canon
    id_data[i,"max_gain_alter"] = max_gain_alter
    id_data[i,"max_loss_alter"] = max_loss_alter
    id_data[i,"max_recip_diff"] = max_recip_diff
  }
}


