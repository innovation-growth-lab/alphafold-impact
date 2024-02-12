"""Project settings. There is no need to edit this file unless you want to change values
from the Kedro defaults. For further information, including these default values, see
https://docs.kedro.org/en/stable/kedro_project_setup/settings.html."""

# Instantiated project hooks.
# For example, after creating a hooks.py and defining a ProjectHooks class there, do
# from pandas_viz.hooks import ProjectHooks

# Hooks are executed in a Last-In-First-Out (LIFO) order.
# HOOKS = (ProjectHooks(),)

# Installed plugins for which to disable hook auto-registration.
# DISABLE_HOOKS_FOR_PLUGINS = ("kedro-viz",)

# Class that manages storing KedroSession data.
from pathlib import Path  # noqa: E402
import itertools

from kedro_viz.integrations.kedro.sqlite_store import SQLiteStore  # noqa: E402

SESSION_STORE_CLASS = SQLiteStore
# Keyword arguments to pass to the `SESSION_STORE_CLASS` constructor.
SESSION_STORE_ARGS = {"path": str(Path(__file__).parents[2])}

# Directory that holds configuration.
# CONF_SOURCE = "conf"

# Class that manages how configuration is loaded.
from kedro.config import OmegaConfigLoader  # noqa: E402

CONFIG_LOADER_CLASS = OmegaConfigLoader
# Keyword arguments to pass to the `CONFIG_LOADER_CLASS` constructor.
CONFIG_LOADER_ARGS = {
    "base_env": "base",
    "default_run_env": "local",
    "config_patterns": {
        "globals": ["parameters*", "parameters*/**", "**/parameters*", "globals*"]
    },
}

DYNAMIC_PIPELINES_MAPPING = {
    "gtr": ["projects", "outcomes/publications", "organisations", "funds"],
    "oa": ["cites", "cited_by"],
    "lens": list(
        itertools.product(
            ["united_states", "european_union"],
            [f"{i:02d}" for i in range(1, 13)],
            [2022, 2023, 2024],
        )
    ),
    "s2": [0,1,2,3,4,5],
    "nsf": ["2018", "2019", "2020", "2021", "2022", "2023", "2024"],
}
