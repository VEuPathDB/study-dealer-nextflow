import java.nio.file.DirectoryStream
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

// Channel.fromPath walks every directory under a glob's literal prefix to the
// pattern's depth, matching or not. This expands one path element at a time and
// only descends into directories whose names match, the way a shell glob does.
class PrunedGlob {

    private static final String GLOB_CHARS = /[*?\[{]/

    static List<Path> expand(String pattern) {
        if(pattern.contains('**')) {
            throw new IllegalArgumentException("PrunedGlob does not support '**': ${pattern}")
        }

        def elements = pattern.split('/') as List
        def firstGlob = elements.findIndexOf { it =~ GLOB_CHARS }
        if(firstGlob < 0) {
            def path = Paths.get(pattern)
            return Files.isRegularFile(path) ? [path] : []
        }

        def prefix = elements.take(firstGlob).join('/')
        def root = Paths.get(prefix ?: (pattern.startsWith('/') ? '/' : '.'))
        def remaining = elements.drop(firstGlob)

        List<Path> current = Files.isDirectory(root) ? [root] : []
        remaining.eachWithIndex { element, i ->
            def isLast = i == remaining.size() - 1
            current = current.collectMany { dir -> children(dir, element, isLast) }
        }

        return current.sort()
    }

    private static List<Path> children(Path dir, String element, boolean wantFiles) {
        if(!(element =~ GLOB_CHARS)) {
            def child = dir.resolve(element)
            return keep(child, wantFiles) ? [child] : []
        }

        def matches = []
        DirectoryStream<Path> stream = Files.newDirectoryStream(dir, element)
        try {
            for(Path child : stream) {
                if(!child.fileName.toString().startsWith('.') && keep(child, wantFiles)) {
                    matches << child
                }
            }
        }
        finally {
            stream.close()
        }
        return matches
    }

    private static boolean keep(Path path, boolean wantFiles) {
        return wantFiles ? Files.isRegularFile(path) : Files.isDirectory(path)
    }
}
