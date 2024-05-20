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
from alphafold_impact.utils.scrape_constants import (
    COUNTRY_CODES,
    US_STATES,
    EXPERIMENTAL_KEYWORDS,
    COMPUTATIONAL_KEYWORDS,
    STRUCTURAL_BIOLOGY_KEYWORDS,
)

logger = logging.getLogger(__name__)


class OpenSyllabusScraper:
    """A class to scrape data from the Open Syllabus website."""

    def __init__(self):
        options = Options()
        options.add_argument("--no-sandbox")
        options.add_argument("--disable-dev-shm-usage")
        prefs = {"download.default_directory": "."}
        options.add_experimental_option("prefs", prefs)

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
    def generate_url(country_name, country_codes, year):
        """Generates the URL to scrape data from the Open Syllabus website."""
        base_url = "https://analytics.opensyllabus.org/record/syllabi"
        params = {
            "field_name": "Biology",
            "field_ids": 9,
            "size": 10000,
            "country_name": country_name,
            "institution_country_codes": country_codes,
            "syllabus_year_start": year,
            "syllabus_year_end": year,
        }
        param_str = "&".join(f"{k}={v}" for k, v in params.items())
        return f"{base_url}?{param_str}"

    @staticmethod
    def generate_state_url(institution_state_codes, state_query, year):
        """Generates the URL to scrape data from the Open Syllabus website."""
        base_url = "https://analytics.opensyllabus.org/record/syllabi"
        params = {
            "field_name": "Biology",
            "field_ids": 9,
            "size": 10000,
            "institution_state_codes": institution_state_codes,
            "state_query": state_query,
            "syllabus_year_start": year,
            "syllabus_year_end": year,
        }
        param_str = "&".join(f"{k}={v}" for k, v in params.items())
        return f"{base_url}?{param_str}"

    @staticmethod
    def generate_keyword_url(keyword):
        """Generates the URL to scrape data from the Open Syllabus website."""
        base_url = "https://analytics.opensyllabus.org/record/syllabi"
        params = {
            "field_name": "Biology",
            "field_ids": 9,
            "size": 10000,
            "syllabus_query": keyword,
        }
        param_str = "&".join(f"{k}={v}" for k, v in params.items())
        return f"{base_url}?{param_str}"

    def collect_data(self, target, country_name, country_codes, year):
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
        logger.info("Collecting data for %s, %s, %s", country_name, country_codes, year)
        if target == "country":
            url = self.generate_url(country_name, country_codes, year)
        else:
            url = self.generate_state_url(country_codes, country_name, year)

        for _ in range(2):  # Retry up to 2 times
            try:
                self.driver.get(url)

                # Wait until the table is loaded
                WebDriverWait(self.driver, 300).until(
                    EC.presence_of_element_located((By.TAG_NAME, "tr"))
                )
                sleep(random.randint(1, 45))

                html = self.driver.page_source
                soup = BeautifulSoup(html, "html.parser")
                data = []
                for row in soup.find_all("tr"):
                    cells = row.find_all("td")
                    if len(cells) > 1:
                        try:
                            course_meta = cells[1].find_all("p")
                            course_code = course_meta[0].text
                            title = course_meta[1].text
                            university_year = course_meta[2].text.rsplit(", ", 1)
                            university = university_year[0]
                            year = int(university_year[1])
                            singleton_href = cells[1].find("a")["href"]
                            data.append(
                                {
                                    "course_code": course_code,
                                    "title": title,
                                    "university": university,
                                    "year": year,
                                    "singleton_href": singleton_href,
                                }
                            )
                        except Exception as e:  # pylint: disable=broad-except
                            logger.error("Error parsing row: %s", e)
                return data
            except Exception as e:  # pylint: disable=broad-except
                logger.error("Error collecting data: %s", e)

        logger.error("Failed to collect data after 2 attempts")
        return []

    def collect_keyword_data(self, keyword):
        """
        Collects data from the Open Syllabus website for the specified
            keyword.

        Args:
            keyword (str): The keyword to search for.

        Returns:
            list: A list of dictionaries containing the collected data.
        """
        logger.info("Collecting data for %s", keyword)

        url = self.generate_keyword_url(keyword)

        for _ in range(2):  # Retry up to 2 times
            try:
                self.driver.get(url)

                # Wait until the table is loaded
                WebDriverWait(self.driver, 300).until(
                    EC.presence_of_element_located((By.TAG_NAME, "tr"))
                )
                sleep(random.randint(1, 45))

                html = self.driver.page_source
                soup = BeautifulSoup(html, "html.parser")
                data = []
                for row in soup.find_all("tr"):
                    cells = row.find_all("td")
                    if len(cells) > 1:
                        try:
                            course_meta = cells[1].find_all("p")
                            course_code = course_meta[0].text
                            title = course_meta[1].text
                            university_year = course_meta[2].text.rsplit(", ", 1)
                            university = university_year[0]
                            year = int(university_year[1])
                            singleton_href = cells[1].find("a")["href"]
                            text_snippet = course_meta[3].text
                            data.append(
                                {
                                    "course_code": course_code,
                                    "title": title,
                                    "university": university,
                                    "year": year,
                                    "singleton_href": singleton_href,
                                    "text_snippet": text_snippet,
                                }
                            )
                        except Exception as e:  # pylint: disable=broad-except
                            logger.error("Error parsing row: %s", e)
                return data
            except Exception as e:  # pylint: disable=broad-except
                logger.error("Error collecting data: %s", e)

        logger.error("Failed to collect data after 2 attempts")
        return []


def main(years, target, output_file_prefix):
    """Scrape data from Open Syllabus for the specified years and country codes."""
    scraper = OpenSyllabusScraper()
    scraper.prepare_page(
        "https://analytics.opensyllabus.org/verify?token=69fd974e-fe50-49f4-9016-6dbd9fccfc64"
    )

    if target == "country":
        for year, country_dict in itertools.product(years, COUNTRY_CODES):
            print(year, country_dict["country_name"])
            country_data = scraper.collect_data(
                target,
                country_dict["country_name"],
                country_dict["country_code"],
                year,
            )
            for item in country_data:
                item["country"] = country_dict["country_name"]
            df = pd.DataFrame(country_data)
            output_file = (
                f"{output_file_prefix}_{country_dict['country_name']}_{year}.parquet"
            )
            df.to_parquet("~/projects/alphafold-impact/data/01_raw/os/" + output_file)
    elif target == "state":
        for year, state_dict in itertools.product(years, US_STATES):
            print(year, state_dict["state_query"])
            state_data = scraper.collect_data(
                target,
                state_dict["state_query"],
                state_dict["institution_state_codes"],
                year,
            )
            for item in state_data:
                item["state"] = state_dict["state_query"]
            df = pd.DataFrame(state_data)
            output_file = (
                f"{output_file_prefix}_{state_dict['state_query']}_{year}.parquet"
            )
            df.to_parquet("~/projects/alphafold-impact/data/01_raw/os/" + output_file)
    elif target == "keyword":
        for keyword in (
            EXPERIMENTAL_KEYWORDS + COMPUTATIONAL_KEYWORDS + STRUCTURAL_BIOLOGY_KEYWORDS
        ):
            print(keyword)
            keyword_data = scraper.collect_keyword_data(keyword)
            df = pd.DataFrame(keyword_data)
            output_file = f"{output_file_prefix}_{keyword}.parquet"
            df.to_parquet(
                "~/projects/alphafold-impact/data/01_raw/os_kw/" + output_file
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scrape data from Open Syllabus.")

    parser.add_argument(
        "--years",
        nargs="+",
        type=int,
        help="List of years to scrape data for.",
    )
    parser.add_argument(
        "--output_file_prefix",
        type=str,
        help="Prefix for output file names.",
        default="syllabii",
    )
    # create argument for choosing either country or US state
    parser.add_argument(
        "--target",
        type=str,
        help="Scrape data for country, state, or keyword.",
    )
    args = parser.parse_args()

    main(args.years, args.target, args.output_file_prefix)
