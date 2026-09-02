#!/bin/bash

# Determine project root relative to the script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Define paths relative to project root
SRC_DIR="$PROJECT_ROOT/.agents/skills"
CLAUDE_DEST="$PROJECT_ROOT/.claude/skills"
CODEX_DEST="$PROJECT_ROOT/.codex/skills"

echo "🔄 Syncing Agent Skills across Claude and Codex environments..."

# Ensure target directories exist
mkdir -p "$CLAUDE_DEST"
mkdir -p "$CODEX_DEST"

# Iterate over files in the source skills directory
if [ -d "$SRC_DIR" ]; then
    for file_path in "$SRC_DIR"/*; do
        # Ensure it's a file
        [ -f "$file_path" ] || continue
        
        filename=$(basename "$file_path")
        
        # Create relative symlinks to prevent broken paths on different machines
        ln -sf "../../.agents/skills/$filename" "$CLAUDE_DEST/$filename"
        ln -sf "../../.agents/skills/$filename" "$CODEX_DEST/$filename"
    done
    echo "✅ Symlinks successfully generated in .claude/skills/ and .codex/skills/"
else
    echo "❌ Error: Source directory $SRC_DIR does not exist."
    exit 1
fi

