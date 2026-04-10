#!/usr/bin/env python3
"""
ChEMBL Drug Information Pipeline
Fetches: Drug ID, Name, Approval Year, Max Phase, 
         Mechanism of Action, Action Type, Disease


- max_phase=4 shows approval year or "Approved (Year Unknown)"
- Searches BOTH mesh_heading AND efo_term for disease
- Approval column shows YEAR only (e.g., "2015", "2003")
"""

"""
Command line usage:
    python chembl_drugs_searching_Sarib.py 
        --gene EGFR 
        --disease "Colorectal Neoplasms" 
        --max-drugs 10 
        --output EFGR_Drugs.xlsx 
        --logs

"""
import time
import argparse
import pandas as pd
from chembl_webresource_client.new_client import new_client


# Setup logging - only shows when --logs flag is used

LOGGING_ENABLED = False

def log_message(msg, level="INFO"):
    """Simple logging function - only prints if LOGGING_ENABLED is True"""
    if LOGGING_ENABLED:
        timestamp = time.strftime("%H:%M:%S")
        print(f"[{timestamp}] [{level}] {msg}")


# Step 1: Find ChEMBL target ID for a given gene

def find_target(gene_symbol):
    """Search for target by gene symbol, return ChEMBL target ID."""
    log_message(f"Searching for target: {gene_symbol}")
    target_api = new_client.target
    
    results = target_api.search(gene_symbol)
    
    if not results:
        log_message(f"No target found for: {gene_symbol}", "WARNING")
        return None
    
    for item in results:
        if item.get("target_type") == "SINGLE PROTEIN" and "Homo sapiens" in item.get("organism", ""):
            chembl_id = item["target_chembl_id"]
            log_message(f"Found human target: {item.get('pref_name', 'Unknown')} ({chembl_id})")
            return chembl_id
    
    chembl_id = results[0]["target_chembl_id"]
    log_message(f"Using fallback target: {results[0].get('pref_name', 'Unknown')} ({chembl_id})", "WARNING")
    return chembl_id


# Step 2: Get mechanisms (drugs) for this target

def get_mechanisms_for_target(target_id, max_drugs=50):
    """Get all mechanisms (drug-target interactions) for a target."""
    log_message(f"Retrieving mechanisms for target: {target_id}")
    mech_api = new_client.mechanism
    
    mechanisms = mech_api.filter(target_chembl_id=target_id)
    mech_list = list(mechanisms[:max_drugs])
    
    log_message(f"Found {len(mech_list)} mechanisms")
    return mech_list


# Step 3: Get approval year for a drug (FIXED LOGIC)

def get_approval_year(mol_data):
    """
    Extract approval year from molecule data.
    If max_phase=4 (approved), returns year or "Approved (Year Unknown)".
    If max_phase < 4, returns "Not Approved".
    """
    max_phase = float(mol_data.get("max_phase", -1))
    
    # Case 1: Drug is NOT approved (max_phase < 4)
    if max_phase < 4:
        return "Not Approved"
    
    # Case 2: Drug IS approved (max_phase == 4)
    if max_phase == 4:
        # Try first_approval field (usually a year as integer)
        first_approval = mol_data.get("first_approval")
        if first_approval:
            # Convert to string, handle if it's a full date
            year_str = str(first_approval)
            # If it's a full date like "2015-06-24", extract just the year
            if "-" in year_str:
                year_str = year_str.split("-")[0]
            log_message(f"    Found approval year: {year_str}", "DEBUG")
            return year_str
        
        # Try usan_year (USAN approval year)
        usan_year = mol_data.get("usan_year")
        if usan_year:
            log_message(f"    Found USAN year: {usan_year}", "DEBUG")
            return str(usan_year)
        
        # Try year_of_approval
        year_of_approval = mol_data.get("year_of_approval")
        if year_of_approval:
            log_message(f"    Found year_of_approval: {year_of_approval}", "DEBUG")
            return str(year_of_approval)
        
        # If no year found but drug is approved
        log_message(f"    No approval year found for approved drug", "WARNING")
        return "Approved (Year Unknown)"
    
    # Fallback for any other max_phase value
    return "Not Approved"


# Step 4: Fetch drug details (Name, Approval Year, Max Phase, Disease)

def get_drug_details(molecule_id, api_delay=0.1):
    """
    Fetch: Name, Approval Year, Max Phase, Disease
    Searches BOTH mesh_heading AND efo_term for disease.
    """
    time.sleep(api_delay)
    
    molecule_api = new_client.molecule
    indication_api = new_client.drug_indication
    
    # Get molecule info
    mol = molecule_api.get(molecule_id)
    if not mol:
        log_message(f"Failed to get molecule: {molecule_id}", "ERROR")
        return None
    
    # Name
    name = mol.get("pref_name", "N/A")
    if name == "N/A":
        synonyms = mol.get("molecule_synonyms", [])
        if synonyms:
            name = synonyms[0].get("molecule_synonym", "N/A")
    
    # Max phase (numeric)
    max_phase = mol.get("max_phase", "N/A")
    
    # Approval Year (FIXED - now handles max_phase=4 correctly)
    approval_year = get_approval_year(mol)
    
    # Disease indications - searches BOTH mesh_heading AND efo_term
    disease_str = "N/A"
    
    try:
        indications = indication_api.filter(molecule_chembl_id=molecule_id)
        ind_list = list(indications)
        log_message(f"  indication_api returned {len(ind_list)} records for {molecule_id}", "DEBUG")
        
        if len(ind_list) > 0:
            diseases = []
            for ind in ind_list[:15]:
                # Search mesh_heading (MeSH terms)
                if ind.get("mesh_heading"):
                    diseases.append(ind["mesh_heading"])
                    log_message(f"    Found mesh_heading: {ind['mesh_heading']}", "DEBUG")
                
                # ALSO search efo_term (Experimental Factor Ontology)
                if ind.get("efo_term"):
                    diseases.append(ind["efo_term"])
                    log_message(f"    Found efo_term: {ind['efo_term']}", "DEBUG")
                
                # Fallback fields
                elif ind.get("disease_mesh_name"):
                    diseases.append(ind["disease_mesh_name"])
                elif ind.get("disease_efo_term"):
                    diseases.append(ind["disease_efo_term"])
            
            if diseases:
                disease_str = "; ".join(set(diseases))
                log_message(f"  ✓ Disease captured for {name}: {len(set(diseases))} unique terms", "DEBUG")
        
    except Exception as e:
        log_message(f"Error fetching disease for {molecule_id}: {str(e)}", "ERROR")
    
    return {
        "Drug_ID": molecule_id,
        "Name": name,
        "Approval_Year": approval_year,  # Now shows year or appropriate status
        "Max_Phase": max_phase,
        "Disease": disease_str
    }


# Main pipeline

def main():
    global LOGGING_ENABLED
    print("Sarib Manzoor - ChEMBL Drug Information Pipeline")
    parser = argparse.ArgumentParser(
        description="Fetch drug information from ChEMBL",
        epilog="Example: python chembl_drugs.py --gene EGFR --disease 'Colorectal Neoplasms' --max-drugs 30"
    )
    
    parser.add_argument("--gene", "-g", type=str, required=True,
                        help="Target gene symbol (e.g., TP53, EGFR)")
    parser.add_argument("--disease", "-d", type=str, default=None,
                        help="Disease name for filtering (optional)")
    parser.add_argument("--output", "-o", type=str, default="chembl_drugs.xlsx",
                        help="Output Excel file name")
    parser.add_argument("--max-drugs", type=int, default=30,
                        help="Maximum number of drugs to fetch (default: 30)")
    parser.add_argument("--delay", type=float, default=0.1,
                        help="Seconds between API calls (default: 0.1)")
    parser.add_argument("--logs", action="store_true",
                        help="Show detailed debug logging output")
    
    args = parser.parse_args()
    
    LOGGING_ENABLED = args.logs
    
    print("=" * 70)
    print("ChEMBL Drug Pipeline")
    print("=" * 70)
    print(f"Target gene: {args.gene}")
    print(f"Max drugs  : {args.max_drugs}")
    print(f"Disease    : {args.disease if args.disease else 'No filter'}")
    print(f"Columns    : Drug_ID, Name, Approval_Year, Max_Phase, Mechanism_of_Action, Action_Type, Disease")
    if args.logs:
        print(f"Logging    : Enabled")
    print("=" * 70)
    
    start_time = time.time()
    
    # Step 1: Find target
    target_id = find_target(args.gene)
    if not target_id:
        print("ERROR: No target found. Exiting.")
        return
    
    # Step 2: Get mechanisms
    mechanisms = get_mechanisms_for_target(target_id, max_drugs=args.max_drugs)
    
    if not mechanisms:
        print("ERROR: No mechanisms found for this target.")
        print("this could be due to the target not having direct drug interactions")
        return
    
    # Step 3: Build results
    print(f"\nFetching drug details...")
    
    results = []
    for idx, mech in enumerate(mechanisms):
        mol_id = mech.get("molecule_chembl_id")
        if not mol_id:
            continue
        
        log_message(f"Processing {idx+1}/{len(mechanisms)}: {mol_id}", "INFO")
        
        # Get mechanism-specific fields
        moa = mech.get("mechanism_of_action", "N/A")
        action_type = mech.get("action_type", "N/A")
        
        # Get drug details
        drug_info = get_drug_details(mol_id, api_delay=args.delay)
        
        if drug_info:
            drug_info["Mechanism_of_Action"] = moa
            drug_info["Action_Type"] = action_type
            results.append(drug_info)
    
    # Optional: Filter by disease
    if args.disease and results:
        original = len(results)
        filtered = []
        for drug in results:
            if drug["Disease"] == "N/A" or args.disease.lower() in drug["Disease"].lower():
                filtered.append(drug)
        results = filtered
        if args.logs:
            print(f"\nDisease filter: {len(results)}/{original} drugs retained")
    
    # Save to Excel
    elapsed_time = time.time() - start_time
    
    if results:
        df = pd.DataFrame(results)
        
        # Updated columns
        columns = ["Drug_ID", "Name", "Approval_Year", "Max_Phase", 
                   "Mechanism_of_Action", "Action_Type", "Disease"]
        df = df[columns]
        
        df.to_excel(args.output, index=False)
        
        print("\n" + "=" * 70)
        print(f"SUCCESS! Completed in {elapsed_time:.1f} seconds")
        print(f"   Saved {len(results)} drugs to {args.output}")
        print("=" * 70)
        
        # Summary
        approved_with_year = sum(1 for d in results if d["Approval_Year"] and d["Approval_Year"].isdigit())
        approved_no_year = sum(1 for d in results if d["Approval_Year"] == "Approved (Year Unknown)")
        not_approved = sum(1 for d in results if d["Approval_Year"] == "Not Approved")
        
        print(f"\nSummary:")
        print(f"  Total drugs: {len(results)}")
        print(f"  Approved with year: {approved_with_year}")
        print(f"  Approved (year unknown): {approved_no_year}")
        print(f"  Not approved / in trials: {not_approved}")
        print(f"  Drugs with disease info: {sum(1 for d in results if d['Disease'] != 'N/A')}/{len(results)}")
    else:
        print("\n⚠ No drugs found matching criteria.")

if __name__ == "__main__":
    main()