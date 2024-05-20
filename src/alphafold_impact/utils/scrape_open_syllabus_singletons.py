# %%
import logging
import random
import itertools
from time import sleep
import argparse
import pandas as pd
from selenium import webdriver
from bs4 import BeautifulSoup
from selenium.webdriver.chrome.service import (  # pylint: disable=ungrouped-imports
    Service,
)
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from webdriver_manager.chrome import ChromeDriverManager


logger = logging.getLogger(__name__)


class OpenSyllabusScraper:
    """A class to scrape data from the Open Syllabus website."""

    def __init__(self):
        options = Options()
        options.add_argument("--no-sandbox")
        options.add_argument("--disable-dev-shm-usage")
        prefs = {"download.default_directory": "."}
        options.add_experimental_option("prefs", prefs)
        options.add_argument("user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3")

        # Get instance of web driver
        self.driver = webdriver.Chrome(
            service=Service(ChromeDriverManager().install()), options=options
        )

    def prepare_page(self, url):
        """
        Opens the specified URL in a web browser and prepares the
            page for scraping.

        Args:
            url (str): The URL of the page to be scraped.
        """
        self.driver.get(url)
        sleep(60)

    @staticmethod
    def generate_url(singleton):
        """Generates the URL to scrape data from the Open Syllabus website."""
        return "https://analytics.opensyllabus.org" + singleton

    def collect_data(self, singleton):
        """
        Collects data from the Open Syllabus website for the specified
            country, country codes, and year.

        Args:
            country_name (str): The name of the country.
            country_codes (str): The country codes.
            year (int): The year for which to collect data.

        Returns:
            list: A list of dictionaries containing the collected data.
        """
        logger.info("Collecting data for %s", singleton)
        url = self.generate_url(singleton)

        for _ in range(2):  # Retry up to 2 times
            try:
                self.driver.get(url)
                WebDriverWait(self.driver, 300).until(
                    EC.presence_of_element_located((By.TAG_NAME, "h2"))
                )

                sleep(random.randint(1, 10))

                html = self.driver.page_source
                soup = BeautifulSoup(html, "html.parser")
                data = {}

                # Find all h2 elements
                h2_elements = soup.find_all("h2")

                # Get the previous sibling of the first h2
                previous_sibling = h2_elements[0].find_previous_sibling()
                data["Introduction"] = previous_sibling.text if previous_sibling else ""

                for h2 in h2_elements:
                    # Get the next sibling of h2
                    sibling = h2.find_next_sibling()

                    # Check if the sibling is a div or ul
                    if sibling.name == "div":
                        # Extract the text from the div
                        data[h2.text] = sibling.text
                    elif sibling.name == "ul":
                        # Extract the text from each li in the ul
                        data[h2.text] = [li.text for li in sibling.find_all("li")]
                # for any values that are lists, ", " join them
                for key, value in data.items():
                    if isinstance(value, list):
                        data[key] = ", ".join(value)
            except Exception as e:  # pylint: disable=broad-except
                logger.error("Error collecting data: %s", e)

            if data:
                return data

        logger.error("Failed to collect data after 2 attempts")
        return []


def main():
    """Scrape data from Open Syllabus for the specified years and country codes."""
    scraper = OpenSyllabusScraper()
    scraper.prepare_page(
        "https://analytics.opensyllabus.org/verify?token=8c397740-5a63-40f8-a737-aae64f67323c"
    )
    data_df = pd.read_parquet(
        "~/projects/alphafold-impact/data/02_intermediate/main_keywords.parquet"
    )

    for i, singleton in enumerate(data_df["singleton_href"].unique()):
        if i < 1150:
            continue
        print(i, singleton)
        data = scraper.collect_data(singleton)
        df = pd.DataFrame([data])
        # extract id from singleton
        singleton_id = singleton.split("id=")[-1]
        df.to_parquet(
            f"~/projects/alphafold-impact/data/02_intermediate/os_kw/{singleton_id}.parquet"
        )
        # if i > 1150:
        #     break


if __name__ == "__main__":
    main()
