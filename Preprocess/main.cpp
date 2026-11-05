// mpcns_pre <setup-file | -> [name=value ...]
#include <cstdio>
#include <exception>
#include <string>

#include "Pipeline.h"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::printf("usage: %s <setup-file | -> [more-setup-files] [name=value ...]\n"
                    "  task=generate|inspect|export|refine|demo_fields  (see Preprocess/INPUT/setup.txt)\n", argv[0]);
        return 2;
    }
    try
    {
        uns::Params p;
        // any number of setup files (later ones win), then name=value overrides in order
        for (int i = 1; i < argc; ++i)
        {
            const std::string a = argv[i];
            if (a == "-") continue;
            if (a.find('=') == std::string::npos) p.load_file(a); else p.set_override(a);
        }
        const std::string task = p.get_string("task", "generate");
        if (task == "generate") return pre::run_generate(p);
        if (task == "inspect") return pre::run_inspect(p);
        if (task == "export") return pre::run_export(p);
        if (task == "refine") return pre::run_refine(p);
        if (task == "demo_fields") return pre::run_demo_fields(p);
        std::fprintf(stderr, "unknown task '%s'\n", task.c_str());
        return 2;
    }
    catch (const std::exception &e)
    {
        std::fprintf(stderr, "mpcns_pre: error: %s\n", e.what());
        return 1;
    }
}
