
translate_species_to_hs = function(species= "Rattus norvegicus",
                                   genes){
  require(msigdbr)
  
  m_df = msigdbr(species = "Rattus norvegicus")
  
  m_df_2=m_df %>% filter(gene_symbol %in% genes)%>%
    distinct(gene_symbol,human_gene_symbol)
  
  one_to_none <- genes[!genes %in% m_df_2$gene_symbol]
  length(one_to_none)
  
  m_df_count= m_df_2 %>% 
    group_by(gene_symbol)%>% count()
  
  one_to_one = m_df_count %>% filter(n==1)%>% pull(gene_symbol)
  length(one_to_one)
  
  ambiguous = m_df_count %>% filter(n!=1)%>% pull(gene_symbol)
  length(ambiguous)
  
  print( paste("one to one:", length(one_to_one)))
  print( paste("one to none:", length(one_to_none)))
  print( paste("ambiguous maps:", length(ambiguous)))
  
  return(list("df"= m_df_2,
              "genes"= list("oto"= one_to_one,
                            "otn"=one_to_none,
                            "amb"= ambiguous)
  )
  )
}