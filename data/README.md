# data/

This is intentionally empty. The repository root is itself the BIDS dataset (ds004725 v1.0.2). MRI content is stored in git-annex and is not committed. Fetch it with `scripts/setup_annex.sh` and then `datalad get <path>`. Derived outputs go in `derivatives/`, which is git-ignored for images and for the row-level phenotype frame.
