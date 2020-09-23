# -*- coding: utf-8 -*-

# This file is part of the magpurify2 package, available at:
# https://github.com/apcamargo/magpurify2
#
# Magpurify2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: antoniop.camargo@gmail.com

import logging
import sys
from argparse import Namespace
import magpurify2


def main(args):
    logger = logging.getLogger("timestamp")
    if not args.fast_mode and not args.taxonomy_database:
        logger.error(
            "The `--taxonomy_database` parameter is required when not using "
            "`--fast_mode`."
        )
        sys.exit(1)
    if args.fast_mode:
        logger.info("Executing MAGpurify2 in fast mode.")
    magpurify2.composition.main(
        Namespace(**magpurify2.cli.default_values["composition"], **vars(args))
    )
    magpurify2.coverage.main(
        Namespace(**magpurify2.cli.default_values["coverage"], **vars(args))
    )
    if not args.fast_mode:
        magpurify2.codon_usage.main(
            Namespace(**magpurify2.cli.default_values["codon_usage"], **vars(args))
        )
        magpurify2.taxonomy.main(
            Namespace(**magpurify2.cli.default_values["taxonomy"], **vars(args))
        )
    # Because both the `end_to_end` and the `filter` modules have the `probability_threshold`
    # and `suffix` parameters, the values are removed from the later to avoid conflict
    magpurify2.cli.default_values["filter"].pop("probability_threshold")
    magpurify2.cli.default_values["filter"].pop("suffix")
    magpurify2.filter.main(
        Namespace(**magpurify2.cli.default_values["filter"], **vars(args))
    )
