import streamlit as st
import pyshorteners as pys
import pyperclip
shortner = pys.Shortener()

st.markdown("<h1 style= 'text-align: center'> URL SHORTENER </h1>", unsafe_allow_html= True)
form = st.form("Name")
url = form.text_input("URL HERE")
s_btn=form.form_submit_button("Submit")
def clicked():
    pyperclip.copy(short_url)
    


if s_btn:
    short_url = shortner.tinyurl.short(url)
    st.markdown("<h3>Your short URL</h3>",unsafe_allow_html=True)
    st.markdown(f"<h6 style='text-align: center'>{short_url}</h6>", unsafe_allow_html= True)
    st.button("Copy", on_click= clicked())


