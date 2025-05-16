import os
import re
import time
import requests
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
import streamlit as st
from datetime import datetime, date

# Set Entrez email
Entrez.email = "swetaraibms@gmail.com"

st.title("Search Keyword and extract PubMed Ids linked to Clinical trial IDs")

def sanitize_keyword(keyword):
    return re.sub(r'[^A-Za-z0-9_]', '_', keyword.strip())

def search_geo_accessions(keyword, retmax=10, retstart=0):
    handle = Entrez.esearch(db="gds", term=keyword, retmax=retmax, retstart=retstart)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", []), int(record.get("Count", 0))

def fetch_geo_accession_details(id_list):
    accession_numbers = set()
    if id_list:
        handle = Entrez.efetch(db="gds", id=id_list, rettype="full", retmode="text")
        data = handle.read()
        handle.close()
        pattern = r"GSE\d{1,10}"
        found_accessions = re.findall(pattern, data)
        accession_numbers.update(found_accessions)
    return list(accession_numbers)

def fetch_all_geo_accessions(keyword, max_results=100):
    all_accessions = set()
    retstart = 0
    total_count = None
    pbar = st.progress(0)
    while True:
        geo_ids, total_count = search_geo_accessions(keyword, retmax=max_results, retstart=retstart)
        if not geo_ids:
            break
        accessions = fetch_geo_accession_details(geo_ids)
        all_accessions.update(accessions)
        retstart += max_results
        if retstart >= total_count:
            break
        pbar.progress(min(retstart / total_count, 1.0))
    pbar.empty()
    return list(all_accessions)

def fetch_pubmed_ids_from_geo(accession_list):
    base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

    def fetch_pubmed_data(accession):
        url = f"{base_url}?acc={accession}&view=full"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        pubmed_ids = []
        for a in soup.find_all('a', href=True):
            if 'pubmed' in a['href']:
                match = re.search(r'pubmed/(\d+)', a['href'])
                if match:
                    pubmed_ids.append(match.group(1))  # Only the number
        return {'accession': accession, 'Pubmed_ID': ', '.join(pubmed_ids)}

    results = []
    progress_bar = st.progress(0)
    status_text = st.empty()

    with ThreadPoolExecutor(max_workers=5) as executor:
        for i, result in enumerate(executor.map(fetch_pubmed_data, accession_list), 1):
            results.append(result)
            progress = i / len(accession_list)
            progress_bar.progress(progress)
            status_text.text(f"📄 Processed {i}/{len(accession_list)}: {result['accession']}")

    progress_bar.empty()
    status_text.empty()
    return results

def fetch_pubmed_html(pubmed_id):
    url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.text
    except:
        return None

def get_publication_date(html):
    soup = BeautifulSoup(html, 'html.parser')
    pub_date_div = soup.find('span', class_='cit')
    if pub_date_div:
        text = pub_date_div.get_text()
        match = re.search(r'\d{4} \w{3} \d{1,2}', text) or re.search(r'\d{4}', text)
        if match:
            try:
                return datetime.strptime(match.group(), "%Y %b %d")
            except:
                try:
                    return datetime.strptime(match.group(), "%Y")
                except:
                    return None
    return None

def search_nct_in_abstract(html):
    soup = BeautifulSoup(html, 'html.parser')

    abstract_divs = soup.find_all('div', class_=lambda x: x and 'abstract-content' in x)
    for div in abstract_divs:
        abstract_text = div.get_text(separator=' ', strip=True)
        match = re.search(r'NCT\d{8}', abstract_text)
        if match:
            return match.group(0)

    trial_reg_divs = soup.find_all('div', class_=lambda x: x and 'trial-registration' in x)
    for div in trial_reg_divs:
        reg_text = div.get_text(separator=' ', strip=True)
        match = re.search(r'NCT\d{8}', reg_text)
        if match:
            return match.group(0)

    full_text = soup.get_text(separator=' ', strip=True)
    match = re.search(r'NCT\d{8}', full_text)
    if match:
        return match.group(0)

    return None

def process_pubmed_ids(metadata):
    pubmed_ids = set()
    for entry in metadata:
        if entry['Pubmed_ID']:
            ids = [pid.strip() for pid in entry['Pubmed_ID'].split(',')]
            pubmed_ids.update(ids)

    pubmed_ids = sorted(list(pubmed_ids))
    results = []

    progress_bar = st.empty()
    status_text = st.empty()  # Single line status

    total = len(pubmed_ids)
    for i, pid in enumerate(pubmed_ids, 1):
        status_text.text(f"🔎 Processing PubMed ID {i}/{total}: {pid}")

        progress_percent = int((i / total) * 100)
        progress_html = f"""
        <div style="border: 1px solid #ccc; border-radius: 5px; width: 100%; height: 20px;">
            <div style="
                width: {progress_percent}%;
                height: 100%;
                background-color: #007bff;
                border-radius: 5px;
                transition: width 0.3s ease-in-out;">
            </div>
        </div>
        <p style="text-align:center; margin: 0;">{progress_percent}%</p>
        """
        progress_bar.markdown(progress_html, unsafe_allow_html=True)

        html = fetch_pubmed_html(pid)
        if html:
            nct = search_nct_in_abstract(html)
            pub_date = get_publication_date(html)
            if isinstance(pub_date, datetime):
                pub_date = pub_date.date()
            results.append({'Pubmed_ID': pid,
                            'NCT Number': nct if nct else 'NCT Not Found',
                            'Publication_Date': pub_date})
        else:
            results.append({'Pubmed_ID': pid,
                            'NCT Number': 'NCT Not Found',
                            'Publication_Date': None})
        time.sleep(0.5)

    progress_bar.empty()
    status_text.empty()

    df = pd.DataFrame(results)
    df["All Pubmed_IDs"] = ", ".join(pubmed_ids)
    return df

def filter_nct(df, start_date, end_date):
    df["Publication_Date"] = pd.to_datetime(df["Publication_Date"], errors='coerce').dt.date
    st.write(f"Filtering between {start_date} and {end_date}")
    st.write(f"Sample publication dates:\n{df['Publication_Date'].dropna().head()}")
    mask = (df["Publication_Date"].notna()) & \
           (df["Publication_Date"] >= start_date) & \
           (df["Publication_Date"] <= end_date)
    filtered_df = df.loc[mask & df['NCT Number'].str.match(r'^NCT\d+$', na=False)]
    st.write(f"Filtered {len(filtered_df)} records with valid NCTs and dates.")
    return filtered_df

# ========== STREAMLIT UI ==========
keyword = st.text_input("Enter keyword for GEO search:")
start_date = st.date_input("Start publication date", key='start_date')
end_date = st.date_input("End publication date", key='end_date')

# After user clicks Run Search button
if st.button("🚀 Run Search"):
    if not keyword:
        st.error("❗ Please enter a keyword to search.")
    elif not start_date or not end_date:
        st.error("❗ Please select both start and end publication dates.")
    else:
        # Store results in session state to persist across reruns
        if 'accessions' not in st.session_state:
            with st.spinner("Fetching GEO accessions..."):
                st.session_state.accessions = fetch_all_geo_accessions(keyword)
        accessions = st.session_state.accessions
        st.success(f"Found {len(accessions)} GEO accessions.")
        geo_df = pd.DataFrame(accessions, columns=["GEO Accession Number"])
        st.dataframe(geo_df)

        if 'metadata' not in st.session_state:
            with st.spinner("Fetching PubMed metadata..."):
                st.session_state.metadata = fetch_pubmed_ids_from_geo(accessions)
        metadata = st.session_state.metadata
        metadata_df = pd.DataFrame(metadata)
        st.dataframe(metadata_df)

        if 'nct_df' not in st.session_state:
            with st.spinner("Extracting NCT Numbers & Publication Dates..."):
                st.session_state.nct_df = process_pubmed_ids(metadata)
        nct_df = st.session_state.nct_df
        st.dataframe(nct_df)

        if 'filtered_nct_df' not in st.session_state:
            with st.spinner("Filtering based on publication date..."):
                st.session_state.filtered_nct_df = filter_nct(nct_df, start_date, end_date)
        filtered_nct_df = st.session_state.filtered_nct_df
        st.success(f"Found {len(filtered_nct_df)} filtered NCT entries.")
        st.dataframe(filtered_nct_df)

        # Download buttons — just use session state data
        today_str = datetime.today().strftime("%Y-%m-%d")

        # Prepare Excel buffer as before
        import io

        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            metadata_df.to_excel(writer, index=False, sheet_name='Metadata')
            nct_df.to_excel(writer, index=False, sheet_name='NCT Extraction')
            filtered_nct_df.to_excel(writer, index=False, sheet_name='Filtered NCT')

        output.seek(0)  

        st.download_button(
            label="Download All Data (Excel)",
            data=output,
            file_name=f"geo_pubmed_all_data_{sanitize_keyword(keyword)}_{today_str}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
