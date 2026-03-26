
#set seed for stability and reproducibility
set.seed(78577)

#libraries
library(tidyverse)

########### Get input files ###############################

#### Download input files
base_dir <- "C:\\data\\bartz_project"
data_dir <- file.path(base_dir,"data")

#If the data directory doesnt exist, create it
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

#store the zip in a temp directory
tempd <- tempdir()

### If the files dont exist, download them and extract them to the data directory
if( !file.exists(file.path(data_dir,"2023.q1.by_size.csv"))){
  download.file("https://data.bls.gov/cew/data/files/2023/csv/2023_q1_by_size.zip",
                destfile=file.path(tempd,"2023.q1.by_size.zip"))
  unzip(zipfile = file.path(tempd,"2023.q1.by_size.zip"),exdir=data_dir)
}

if( !file.exists(file.path(data_dir,"2024.q1.by_size.csv")) ){
  download.file("https://data.bls.gov/cew/data/files/2024/csv/2024_q1_by_size.zip",
                destfile=file.path(tempd,"2024.q1.by_size.zip"))
  unzip(zipfile = file.path(tempd,"2024.q1.by_size.zip"),exdir=data_dir)
}

##########################################

# Generate Frame Using the Files Above

#I only keep the following columns when reading in
keepcols <- c("area_fips", "industry_code", "disclosure_code","agglvl_code",
              "size_code", "year", "qtr", "qtrly_estabs_count",  "month2_emplvl",
              "month3_emplvl")

#get the 2023 dataset
df1 <- read_csv( file.path(data_dir,"2023.q1.by_size.csv") ,
                col_select = all_of(keepcols))|> 
  mutate( state_fips = str_sub(area_fips,1,2) ) |> #Get the state code from area_fips
  filter(agglvl_code == 62, #I only keep the 2 largest industries under total employment 
        state_fips < 72,  #Remove US territories
        is.na(disclosure_code) ) |> #there are some cells that cannot be published, so I get rid of them
  select(-disclosure_code,-year,-qtr,-area_fips,-agglvl_code) |> #
  rename( py_month2_emplvl = month2_emplvl,
          py_month3_emplvl = month3_emplvl,
          py_qtrly_estabs_count = qtrly_estabs_count)

#get the 2024 dataset
df2 <- read_csv(file.path(data_dir,"2024.q1.by_size.csv"),
                col_select = all_of(keepcols) ) |> 
  mutate( state_fips = str_sub(area_fips,1,2) ) |> 
  filter(agglvl_code == 62,
  state_fips < 72,
   is.na(disclosure_code) ) |> 
  select(-disclosure_code,-area_fips,-agglvl_code)


## columns I use to merge 2023 data onto 2024 data
mergeon <- c("state_fips", "industry_code", "size_code")

df <- df2 |> 
  left_join(df1,by=mergeon) |> ##merge on last years data
  drop_na() |> ## I drop the rows taht didnt merge on cleanly due to difference in diclosure from year to year, there were like 4 of them
mutate( ces_size = recode_values(size_code,   # I recode from QCEW to CES size classes                         
        c(1,2)~1,                             # QCEW https://www.bls.gov/cew/classifications/size/size-titles.htm
        3~2,                                  # CES: Exhibit 1, https://www.bls.gov/opub/hom/ces/design.htm
        4~3,
        5~4,
        6~5,
        7~6,
        8~7,
        9~8)) |>
  select(-size_code) |>    #remove QCEW Size code
  group_by(  state_fips, industry_code, ces_size,  year,   qtr) |> #resummarize variables 
  summarise(across(everything(), sum, na.rm = TRUE)) |> ## Resummarize all columns not in group_by
  ungroup() |>  #ungroup to avoid complications later
  mutate( beta_m2m3    = month3_emplvl/month2_emplvl,         #change from month2 to month3 current year
          py_beta_m2m3 = py_month3_emplvl/py_month2_emplvl,   #change from month2 to month3 previous year
          m2_emp  = month2_emplvl/qtrly_estabs_count,         #current year mean for month 2
          m3_emp  = month3_emplvl/qtrly_estabs_count,         #current year mean for month 3
          py_m2_emp = py_month2_emplvl/py_qtrly_estabs_count, #previous year mean for month 2
          py_m3_emp = py_month3_emplvl/py_qtrly_estabs_count, #previous year mean for month 3
          ces_size = as_factor(ces_size),              #set as categorical
          beta_yoym2 = month2_emplvl/py_month2_emplvl, #Year over year change for month 2
          beta_yoym3 = month3_emplvl/py_month3_emplvl, #Year over year change for month 3      
          estabs_count_sim = ceiling( qtrly_estabs_count/10), #I divide by 10 cause I dont want to generate 10M units,
          estabs_count_sim = if_else(estabs_count_sim<2,2,estabs_count_sim) ## Make at least 2 units in a stratum
   )

### Poisson has too many 0's so I protect against it
live <- function(x,p){
  if_else(x==0, rbinom(n=length(x),size=1,p=p),x )
}

### function below will iterate through the rows of df and generate a population for each row
## I reindex the months to 1 and 2 for the sake of our simulated population

generate <- function( state_fips, industry_code, ces_size, year, qtr, py_qtrly_estabs_count,  qtrly_estabs_count,                
                      month2_emplvl, month3_emplvl, 
                      py_month2_emplvl, py_month3_emplvl, beta_m2m3,  
                      py_beta_m2m3,  m2_emp, m3_emp,                
                      py_m2_emp, py_m3_emp, beta_yoym2, beta_yoym3, estabs_count_sim  ){
     
    p=.95
    #Generate the population
    m1py = rpois(n=estabs_count_sim,lambda=py_m2_emp)
    m1py = live(m1py,p)
  
    m2py = sapply(m1py,\(x){rpois(n=1,lambda=x*py_beta_m2m3)})
    m2py = live(m2py,p)
    
    m1y  = sapply(m1py,\(x){rpois(n=1,lambda=x*beta_yoym2)})
    m1y  = live(m1y,p)
  
    m2y  = sapply(m1y,\(x){rpois(n=1,lambda=x*beta_m2m3)})
    m2y  = live(m2y,p)

    #return a dataframe to the user
    return( tibble(state_fips=state_fips, 
                   industry_code=industry_code,
                   ces_size = ces_size,
                   m1py=m1py,
                   m2py=m2py,
                   m1y = m1y,
                   m2y = m2y)
    )

}

### for each row of the data frame we created, generate samples from it
frame <- df |> pmap(generate) |> bind_rows()


#frame |> summarize(m1y=sum(m1y),m2y=sum(m2y))
#######################################################################################

## Create Sampling Design
# I do a stratified simple random sample and use Neyman Allcation to allocate the sample

#Frame Size
N <- nrow(frame)

#Total Sample Size. I use roughly the proportion of sample CES collects in terms of establishments
n = floor(0.075*N)

### We will need to reiterate this a few times when stratum sample > stratum count

## calculate allocation, this is a temporary holder because 
## we will need to iterate to allocate sample
allocation_temp <- 
frame |> 
  group_by(state_fips, industry_code, ces_size) |> 
  summarize(N_h = n(),
            R = sum(m2py)/sum(m1py),
            Sm1 = sd(m1py),
            Sm2 = sd(m2py),
            rho = cor(m1py,m2py),
            S = Sm1^2 + R^2*Sm2^2 - 2*R*Sm1*Sm2
) |> ungroup()

#### iterative reallocation until we assign all sample

## new_n is used to keep track of the amount of sample left
## to allocate
new_n <- n

#bucket to store allocated sample during loop
allocation <- tibble()

#while we have not assigned all sample run loop
while(new_n > 0){
  ### allocate according to neyman formula
  allocation_temp <- 
  allocation_temp  |> 
    mutate(n_h = ceiling(new_n*N_h*S/(sum(N_h*S))) )

  #check if any n_h is greater than stratum size
  if(any(allocation_temp$n_h > allocation_temp$N_h)){
    
    #if it is, set the sample size to stratum size
    temp <- 
    allocation_temp |> filter(n_h>N_h) |> mutate(n_h=N_h)
    
   #add these items to the final allocation table
    allocation <- temp |> bind_rows(allocation)

    #remove this allocated sample
    new_n <- new_n - sum(temp$n_h)

    #remove the strata we have allocated sample to and repeat if necessary
    allocation_temp <- allocation_temp |> filter(n_h<=N_h)
  } else {
    
    #if we didnt over allocate, then we are done

    #remove the remaining sample we allocated
    new_n <- new_n - sum(allocation_temp$n_h)
    
    #add these strata to the final allocation table
    allocation <- allocation_temp |> bind_rows(allocation)

  }
}

### Create the weights and sort
allocation <- allocation |> 
              mutate(w=N_h/n_h) |> 
              arrange(state_fips,industry_code,ces_size) |> 
  select(state_fips,industry_code,ces_size,n_h,w)

############################################
###### Begin sampling 

#I nest the data frame and left join sample sizes onto it
frame_nest <- 
frame |> 
  group_by(state_fips,industry_code,ces_size) |> 
  nest() |> 
  left_join(allocation,by=c("state_fips", "industry_code", "ces_size"))

##number of samples to dr
n_samps <- 150

#bucket to hold sample estimates
#I will add to this each time I draw a sample and estimate
sample_ests <- tibble()

#draw n_samps and create an estimate for each
for(i in 1:n_samps){
  print(i)

  #run a sample
  s1<- 
  frame_nest |> 
    mutate(sample = map2(data,n_h, \(x,y){slice_sample(x,n=y)})) |> 
    select(-data) |> 
    unnest(sample) |> 
    ungroup()

  #create a horvitz thompsons estimate for a total for month1 and month2
  #In the last step, I add this as a row onto our dataframe sample_ests
  sample_ests <- s1 |> summarize(m1y_est=sum(w*m1y),
                     m2y_est=sum(w*m2y)) |> bind_rows(sample_ests)
}

#I grab the true frame total for each month cause I want to plot it
pop_totals <- frame |> summarize(m1y = sum(m1y), m2y=sum(m2y))

#make scatter plot of estimates and add the population value
sample_ests |> ggplot(aes(x=m1y_est,y=m2y_est)) + 
                geom_point() + 
                geom_point(aes(x=m1y,y=m2y),size=8,color="red",data=pop_totals)

