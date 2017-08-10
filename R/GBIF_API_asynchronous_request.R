



function gbifapi { 

curl -i –user popple_1500:Popple1500 -H "Content-Type: 
application/json" -H "Accept: application/json" -X POST -d "{\"creator\":
\”yourGbifUserName\", \"notification_address\": [\”hugh.burley@mq.edu.au\"], 
\"predicate\": {\"type\":\"and\", \"predicates\": [{\"type\":\"equals\",\"key\":
\"HAS_COORDINATE\",\"value\":\"true\"}, {\"type\":\"equals\", \"key\":
\"TAX O N _KEY\", \"value\":\"$1\"}] }}" 
http://api.gbif.org/v1/occurrence/download/request  >> log.txt 
echo -e "\r\n$1 $2\r\n\r\n----------------\r\n\r\n" >> log.txt
 
}  

$ gbifapi 4140730 "Betula pendula" 
$ gbifapi 2882316 "Fagus sylvatica" 
$ gbifapi 3172358 "Fraxinus excelsior" 
$ gbifapi 5289784 "Hedera helix" 
$ gbifapi 2878688 "Quercus robur"  

# Downloading species occurrence data using the GBIF web-service API (PDF Download Available). 
# Available from: https://www.researchgate.net/publication/308802095_Downloading_species_occurrence_data_using_the_GBIF_web-service_API [accessed Aug 10, 2017].



