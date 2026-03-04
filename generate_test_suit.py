

import pandas as pd
import os
import pickle
import json
from collections import defaultdict

def load_pkl_file(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def save_pkl_file(data, file_path):
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)
    print(f"Updated data saved to {file_path}")
    print(f"Data Count: {len(data)}")
    print( "-----------------" )

def read_ipr_mapping( loc ):

    ipr_map_info_loc = loc + "ipr_mapping.json"

    with open(ipr_map_info_loc, "r", encoding="utf-8") as f:
        ipr_mapping = json.load(f)

    return ipr_mapping

def view_seq_dic( sequence_sets ):

    # print(sequence_sets)

    for uniprot_id, data in sequence_sets.items():
        iprs = data.get("interpro", [])
        sequence = data.get("sequence", "")
        ipr_info_intersection_type = data.get("ipr_info_intersection_type", "")
        iprs_not_in_train = data.get("iprs_not_in_train", "")
    
        print(f"IPRS:\n{iprs}")
        print(f"ipr_info_intersection_type:\n{ipr_info_intersection_type}")
        print(f"iprs_not_in_train:\n{iprs_not_in_train}")
        print(f"Sequence:\n{sequence}")
        print("--------------")
    

def get_sequence( new_dataset ):
    sequences = {}
     
    for entry in new_dataset:
        entry_id = entry.get('uniprot_id', None)
        ipr_numbers = entry.get('ipr_numbers', [])
        sequence = entry.get('sequence', None)

        sequences[entry_id] = {
                "sequence": sequence,
                "length": len(sequence),
                "interpro": ipr_numbers
            }
    
    return sequences

def merge_dic(d1, d2):
    d = defaultdict(dict)
    for other in [d1, d2]:
        for k, v in other.items():
            d[k] = v
    
    return d

def filter_by_length( sequences, n = 1000 ):
    new_sequences = {key:val for key, val in sequences.items() if (sequences[key]["length"] >= 1000 and sequences[key]["length"] <= 2000)}
    return new_sequences


def filter_by_two_set( sequences, trained_sequence ):

   new_sequences = {key:val for key, val in sequences.items() if key not in list(trained_sequence.keys())}
   return new_sequences


def remove_train_match( sequences, loc = "../../paper/general_dataset_go_ipr/cfpgen_general_dataset/" ):
    train_pkl_loc = loc + "train.pkl"

    dataset = load_pkl_file( train_pkl_loc )
    trained_sequence = get_sequence( dataset )
    trained_sequence = filter_by_length( trained_sequence )

    new_sequences = filter_by_two_set( sequences, trained_sequence )

    return new_sequences


def parse_fasta_to_dict(fasta_loc):
    sequences = {}

    with open(fasta_loc, "r") as fp:
        entry_id = None
        sequence_chunks = []

        for line in fp:
            line = line.strip()

            if line.startswith(">"):
                if entry_id:
                    sequences[entry_id] = {
                        "sequence": "".join(sequence_chunks),
                        "length": None,
                        "interpro": []
                    }

                entry_id = line.split("|")[1]
                sequence_chunks = []
            else:
                sequence_chunks.append(line)

        if entry_id:
            sequences[entry_id] = {
                        "sequence": "".join(sequence_chunks),
                        "length": None,
                        "interpro": []
                    }

    return sequences

def gather_add_info( tsv_loc, sequences ):

    uni_kb_tsv = pd.read_csv( tsv_loc, sep='\t', header=0 )

    for _, row in uni_kb_tsv.iterrows():
        entry = row["Entry"]
        length = row["Length"]
        ipr_list = row["InterPro"].split(";") if pd.notna(row["InterPro"]) else []
        # print(ipr_list)

        if entry in sequences and len(ipr_list) != 0:
            
            sequences[entry]["length"] = length
            sequences[entry]["interpro"] = ipr_list[:-1]
        
        else:
            del sequences[entry]

    return sequences


def return_new_sequence_set( loc = None, file_name = "uniprotkb_AND_reviewed_true_AND_model_o_2026_02_18" ):

    fasta_loc = loc + file_name + ".fasta"
    tsv_loc = loc + file_name + ".tsv"

    sequences = parse_fasta_to_dict(fasta_loc)
    sequences = gather_add_info( tsv_loc, sequences )

    ## debug: view entry from dict
    # view_seq_dic( sequences )
    

    print( f"Sequence Length Before Filter By Length: {len(sequences)}" ) # 3174    
    sequences = filter_by_length(sequences)
    print( f"Sequence Length after Filter  By Length: {len(sequences)}" ) # 2093

    sequences = remove_train_match( sequences=sequences )
    print( f"Sequence Length After train match: {len(sequences)}" ) #2092

    return sequences

def filter_cfp_test( loc ):
    test_loc = loc + "test.pkl"

    dataset = load_pkl_file( test_loc )
    sequences = get_sequence( dataset )

    print( f"Sequence Length Before Filter By Length: {len(sequences)}" ) # 8309
    sequences = filter_by_length(sequences)
    print( f"Sequence Length after Filter  By Length: {len(sequences)}" ) # 1

    return sequences

def ipr_info_map_match( sequence_sets, ipr_info, name = None ):

    affected_seq = {}

    ## check new_ipr_info_count
    seq_ipr_not_fnd_in_info_map_cnt = 0
    cnt_sq_affect = 0

    unique_set_ipr = set()
    ipr_info_map_lst = list(ipr_info.keys() )

    ## adding all partial ipr-found sequence
    for uniprot_id, data in sequence_sets.items():
        iprs_in_seq = data.get("interpro", [])
        is_seq_affect = False
        for ipr_in_seq in iprs_in_seq:
            unique_set_ipr.add( ipr_in_seq )
            if( ipr_in_seq not in ipr_info_map_lst ):
                is_seq_affect = True
                affected_seq[uniprot_id] = sequence_sets[uniprot_id].copy()
                break
                
        if( is_seq_affect == True ):
            cnt_sq_affect += 1

    
    # affected_seq[uniprot_id_to_remove] = sequence_sets.pop(uniprot_id_to_remove)
    
    
    ## filtering out all absolutely-no-ipr found ones
    # main_list: ipr_info_map_lst ["a", "b", "c"]
    # list to check: iprs_in_seq ["a,", "b"]
    
    no_ipr_fnd_cnt = 0
    uniprot_id_to_remove = []
    for uniprot_id, data in affected_seq.items():
        affected_seq[uniprot_id]["ipr_info_intersection_type"] = None
        iprs_in_seq = data.get("interpro", [])
        if any( e in ipr_info_map_lst for e in iprs_in_seq ):
            affected_seq[uniprot_id]["ipr_info_intersection_type"] = "partial"
            if all( e in ipr_info_map_lst for e in iprs_in_seq ):
                affected_seq[uniprot_id]["ipr_info_intersection_type"] = "full"          
        else: ## if not a single ipr found, then breaks
            uniprot_id_to_remove.append( uniprot_id )
            no_ipr_fnd_cnt += 1

    
    for uniprot_id in uniprot_id_to_remove:
        del affected_seq[uniprot_id]

                
    # for unq_ipr in list(unique_set_ipr):
    #     if( unq_ipr not in list(ipr_info.keys() ) ):
    #        seq_ipr_not_fnd_in_info_map_cnt += 1
    
    
    # print( f"-----Merged" )
    # print( f"Total instance: {len(sequence_sets)}" ) #2093
    # print( f"Total affected instance: {cnt_sq_affect}" ) #2090
    # print( f"Absolutely No-IPR #: {no_ipr_fnd_cnt}" ) #1692
    # print( f"AtLeast One IPR Containing Set Count: {len(affected_seq)}" ) #398/299 (for 1k-2k length)
    # print( f"Total Unique IPR in Train Set: {len(unique_set_ipr)}" )
    # print( f"IPR Beyond train set: {seq_ipr_not_fnd_in_info_map_cnt}" )
    
    return affected_seq



def retrieve_set( sequence_sets, ipr_info ):

    # recieve only one, with updated param
    partial_seq = ipr_info_map_match( sequence_sets, ipr_info )

    # view_seq_dic( partial_seq )

    return partial_seq
    


def align_sets( sequence_sets, ipr_info ):

    """
    Updates Dictionary Key to CFP-GEN expected .pickle type
    
    :param sequence_sets: dictonary
    """

    result = []

    for uniprot_id, data in sequence_sets.items():

        seq_affect = False
        iprs = data.get("interpro", [])

        ipr_info_list = list(ipr_info.keys())

        # only consider ipr found in ipr_info
        new_iprs = [ipr_info[ipr] for ipr in iprs if ipr in ipr_info_list]

        if( seq_affect == False ):

            entry = {
                "uniprot_id": uniprot_id,
                "ipr_numbers": data.get("interpro", []),
                "go_numbers": {
                    "C": [],
                    "F": [],
                    "P": []
                },
                "ipr_mapped": new_iprs,
                "go_f_mapped": [],
                "iprs_in_train": len(new_iprs),
                "iprs_not_in_train": len(iprs) - len(new_iprs),
                "length": data.get("length", None),
                "sequence": data.get("sequence", "")
            }

            result.append(entry)

    # print(result)

    return result




if "__main__":

    ipr_info = read_ipr_mapping( loc = "../dataset/uniprotkb/" )

    # sequences = return_new_sequence_set( loc = "../dataset/uniprotkb/", file_name="uniprotkb_AND_reviewed_true_AND_length_2026_02_18" ) 
    sequences = return_new_sequence_set( loc = "../dataset/uniprotkb/" )    

    ## cfp-gen test (n > 1000)
    dataset_test_sequences = filter_cfp_test( loc = "../../dataset/general_dataset_go_ipr/cfpgen_general_dataset/" )
    
    # ## combination: join two dic
    merged_sequences = merge_dic( sequences, dataset_test_sequences )
    print(f"Merged Sequence Length: {len(merged_sequences)}")

    # ## inference_test
    pure_new_sequences = filter_by_two_set( sequences, dataset_test_sequences )
    print(f"Pure Sequence Length: {len(pure_new_sequences)}")

    # ### Debugging
    # ipr_info_map_match( dataset_test_sequences, ipr_info, "dataset_test_sequences" )
    ipr_info_map_match( merged_sequences, ipr_info, "merged_sequences" )
    # ipr_info_map_match( pure_new_sequences, ipr_info, "pure_new_sequences" )

    ### Create for set
    # 1. Sequence consisting IPR found in info-map (Full) (sz: 2; with test_set)
    # 2. Sequence consisting IPR found in info-map (Partial) (sz: 2092)
    ## Note: I am Not testing with sequences where there is absolutely no IPR found on map_info as I do not know their mapping.

    ## retrieving new match
    full_ipr_containing_set = retrieve_set( merged_sequences, ipr_info )


    # dataset_test_sequences = align_sets( dataset_test_sequences  )
    # merged_sequences = align_sets( merged_sequences  )
    # pure_new_sequences = align_sets( pure_new_sequences  )
    experiment_set = align_sets( full_ipr_containing_set, ipr_info)

    # ## save each to pickle
    # save_pkl_file( dataset_test_sequences, "../../paper/code/ml4lp/data-bin/uniprotKB/cfpgen_general_dataset/"  + "test_only_large.pkl" )
    # save_pkl_file( merged_sequences, "../../paper/code/ml4lp/data-bin/uniprotKB/cfpgen_general_dataset/"  + "test_w_new_large.pkl" )
    # save_pkl_file( pure_new_sequences, "../../paper/code/ml4lp/data-bin/uniprotKB/cfpgen_general_dataset/"  + "only_new_large.pkl")
    save_pkl_file( experiment_set, "../../paper/code/ml4lp/data-bin/uniprotKB/cfpgen_general_dataset/"  + "experiment.pkl")



    