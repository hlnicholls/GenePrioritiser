#!/usr/bin/env bash
set -euo pipefail
# Create an empty directory tree matching the repository 'example/' folder,
# but rooted at the user's configured PROJECT_DIR (from env, config.sh or config.py).
# This lets users change PROJECT_DIR (for example to '/results') and ensure the
# required subfolders (data_preprocessing, machine_learning, etc.) exist before
# running Step1.

usage() {
    cat <<EOF
Usage: $0 [TARGET_PROJECT_DIR]

Determine the project root in this order:
  1) argument TARGET_PROJECT_DIR
  2) environment variable 'projectDir'
  3) the exported PROJECT_DIR from config/config.sh (if present)
  4) the PROJECT_DIR value in config/config.py if explicitly set
  5) default '.' (repository root)

This script will walk the repository 'example/' directory and create the
corresponding directories under the chosen PROJECT_DIR (the leading 'example'
component is stripped). Files are not copied.
EOF
}

# allow override by first argument
if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

TARGET_ARG="${1:-}"

# 1) prefer CLI arg
if [ -n "$TARGET_ARG" ]; then
    PROJECT_DIR="$TARGET_ARG"
else
    # 2) environment variable projectDir
    if [ -n "${projectDir:-}" ]; then
        PROJECT_DIR="$projectDir"
    else
        # 3) try sourcing config/config.sh if present
        if [ -f "$(dirname "$0")/../../config/config.sh" ]; then
            # shellcheck source=/dev/null
            source "$(dirname "$0")/../../config/config.sh"
            PROJECT_DIR="${PROJECT_DIR:-.}"
        else
            # 4) try parsing config/config.py for a literal PROJECT_DIR assignment
            if [ -f "$(dirname "$0")/../../config/config.py" ]; then
                PROJ_LINE=$(grep -E "^PROJECT_DIR\s*=\s*['\"]" "$(dirname "$0")/../../config/config.py" || true)
                if [ -n "$PROJ_LINE" ]; then
                    PROJECT_DIR=$(echo "$PROJ_LINE" | sed -E "s/PROJECT_DIR\s*=\s*['\"](.*)['\"]/\1/")
                else
                    PROJECT_DIR='.'
                fi
            else
                PROJECT_DIR='.'
            fi
        fi
    fi
fi

# normalize to absolute path
PROJECT_DIR_ABS=$(python3 - <<PY
import os,sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PY
"$PROJECT_DIR")

echo "Creating directory skeleton under project root: $PROJECT_DIR_ABS"
mkdir -p "$PROJECT_DIR_ABS"

REPO_ROOT=$(cd "$(dirname "$0")/../.." && pwd)
EXAMPLE_DIR="$REPO_ROOT/example"

if [ ! -d "$EXAMPLE_DIR" ]; then
    echo "No 'example' directory found at $EXAMPLE_DIR. Nothing to create." >&2
    exit 1
fi

count=0
# Walk example tree and create matching directories under PROJECT_DIR_ABS,
# stripping the leading 'example/' component.
while IFS= read -r -d $'\0' dir; do
    rel=${dir#"$EXAMPLE_DIR"}
    # remove leading slash if present
    rel=${rel#/}
    # strip leading 'example/' component if user included it
    rel=${rel#example/}
    target_dir="$PROJECT_DIR_ABS/$rel"
    mkdir -p "$target_dir"
    count=$((count+1))
done < <(find "$EXAMPLE_DIR" -type d -print0)

echo "Created $count directories under $PROJECT_DIR_ABS (mirroring 'example/')."
echo "If you want a different root, re-run with the target path as the first argument."

exit 0
