// mpcns_pre <setup-file | -> [name=value ...]
#include <cstdio>
#include <filesystem>
#include <vector>
#include <exception>
#include <string>

#include "Pipeline.h"

int main(int argc, char **argv)
{
    if (argc < 2 && !std::filesystem::exists("config.yml"))
    {
        std::printf("usage: %s                      (reads ./config.yml)\n"
                    "       %s <config.yml | setup.txt | -> [more-files] [name=value ...]\n"
                    "  task=generate|inspect|export|refine|demo_fields  (see case_template/*/01_Meshgen)\n", argv[0], argv[0]);
        return 2;
    }
    try
    {
        uns::Params p;
        // config.yml (default) and/or setup files (later ones win), then name=value overrides in order
        p.load_args(std::vector<std::string>(argv + 1, argv + argc));
        if (p.has("cost_profile")) // load-balance profile written by the solver's previous run (optional)
        {
            const std::string cp = p.get_string("cost_profile", "");
            if (std::filesystem::exists(cp)) p.load_file(cp);
            else std::printf("note: cost_profile '%s' not found, using the default cost model\n", cp.c_str());
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
