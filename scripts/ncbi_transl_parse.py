#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

import re
import requests
from bs4 import BeautifulSoup


def get_html():
    url = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi'
    return requests.get(url).text


def parse():
    html = get_html()
    parts = re.split(r'<hr\s*/?\s*>', html)[1:]
    for part in parts:
        lines = []
        soup = BeautifulSoup(part, 'lxml')
        h2 = soup.select_one('h2')
        pre = soup.select_one('pre')
        if not h2 or not pre:
            continue
        lines.append('# ' + h2.text.strip())
        lines.extend(x.strip() for x in pre.text.splitlines() if '=' in x)
        lines.append('\n==SEPARATOR==\n')
        yield lines


def parse_and_print():
    print('# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi\n')
    for lines in parse():
        for line in lines:
            print(line)


if __name__ == '__main__':
    parse_and_print()
